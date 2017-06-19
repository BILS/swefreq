import applicationTemplate
import email.mime.multipart
from email.mime.text import MIMEText
import json
import logging
import peewee
import pymongo
import smtplib
import tornado.template as template
import tornado.web

import db
import handlers
import settings


class Query(handlers.UnsafeHandler):
    def make_error_response(self):
        ret_str = ""

        checks = {
                'dataset': lambda x: "" if x == 'SweGen' else "dataset has to be SweGen\n",
                'ref': lambda x: "" if x == 'hg19' else "ref has to be hg19\n",
                'pos': lambda x: "" if x.isdigit() else "pos has to be digit\n",
        }

        for arg in ['chrom', 'pos', 'dataset', 'referenceAllele', 'allele', 'ref']:
            try:
                val = self.get_argument(arg)
                if arg in checks:
                    ret_str += checks[arg](val)
            except:
                ret_str += arg + " is missing\n"
                if arg in checks:
                    ret_str += checks[arg]("")

        dataset = self.get_argument('dataset', 'MISSING')

        return ret_str

    def get(self, *args, **kwargs):
        the_errors = self.make_error_response()
        if len(the_errors) > 0:
            self.set_status(400);
            self.set_header('Content-Type', 'text/plain');
            self.write(the_errors);
            return

        sChr      = self.get_argument('chrom', '').upper()
        iPos      = self.get_argument('pos', '')
        dataset   = self.get_argument('dataset', '')
        referenceAllele = self.get_argument('referenceAllele', '').upper()
        allele = self.get_argument('allele', '').upper()
        reference = self.get_argument('ref', '').upper()

        exists = lookupAllele(sChr, int(iPos), referenceAllele, allele, reference, dataset)

        if self.get_argument('format', '') == 'text':
            self.set_header('Content-Type', 'text/plain')
            self.write(str(exists))
        else:
            self.write({
                'response': {
                    'exists': exists,
                    'observed': 0,
                    'externalUrl': "%s://%s" % ('https', self.request.host),
                    },
                'query': {
                    'chromosome': sChr,
                    'position': iPos,
                    'referenceAllele': referenceAllele,
                    'allele': allele,
                    'dataset': dataset,
                    'reference': reference
                    },
                'beacon': 'swefreq-beacon'
                })

class Info(handlers.UnsafeHandler):
    def get(self, *args, **kwargs):
        query_uri = "%s://%s/query?" % ('https', self.request.host)
        self.write({
            'id': 'swefreq-beacon',
            'name': 'Swefreq Beacon',
            'organization': 'SciLifeLab',
            'api': '0.3',
            #'description': u'Swefreq beacon from NBIS',
            'datasets': [
                {
                    'id': 'SweGen',
                    # 'description': 'Description',
                    # 'size': { 'variants': 1234, 'samples': 12 },
                    # 'data_uses': [] # Data use limitations
                    'reference': 'hg19'
                },
            ],
            'homepage':  "%s://%s" % ('https', self.request.host),
            #'email': u'swefreq-beacon@nbis.se',
            #'auth': 'None', # u'oauth2'
            'queries': [
                query_uri + 'dataset=SweGen&ref=hg19&chrom=1&pos=55500975&referenceAllele=C&allele=T',
                query_uri + 'dataset=SweGen&ref=hg19&chrom=1&pos=55505551&referenceAllele=A&allele=ACTG&format=text',
                query_uri + 'dataset=SweGen&ref=hg19&chrom=2&pos=41936&referenceAllele=AG&allele=A'
                ] #
            })

def lookupAllele(chrom, pos, referenceAllele, allele, reference, dataset):
    """CHeck if an allele is present in the database
    Args:
        chrom: The chromosome, format matches [1-22XY]
        pos: Coordinate within a chromosome. Position is a number and is 0-based
        allele: Any string of nucleotides A,C,T,G
        alternate: Any string of nucleotides A,C,T,G
        reference: The human reference build that was used (currently unused)
        dataset: Dataset to look in (currently used to select Mongo database)
    Returns:
        The string 'true' if the allele was found, otherwise the string 'false'
    """
    client = pymongo.MongoClient(host=settings.mongo_host, port=settings.mongo_port)

    # The name of the dataset in the database is exac as required by the
    # exac browser we are using.
    if dataset == 'SweGen':
        dataset = 'exac'

    mdb = client[dataset]
    mdb.authenticate(settings.mongo_user, settings.mongo_password)

    # Beacon is 0-based, our database is 1-based in coords.
    pos += 1
    res = mdb.variants.find({'chrom': chrom, 'pos': pos})
    for r in res:
        if r['alt'] == allele and r['ref'] == referenceAllele:
            return True

    return False

class Home(handlers.UnsafeHandler):
    def get(self, *args, **kwargs):
        t = template.Template(applicationTemplate.index)

        has_access = self.is_authorized()
        is_admin   = self.is_admin()

        name = None
        email = None
        if self.current_user:
            name = self.current_user.name
            email = self.current_user.email

        self.write(t.generate(user_name  = name,
                              has_access = has_access,
                              email      = email,
                              is_admin   = is_admin,
                              ExAC       = settings.exac_server))


class GetDataset(handlers.UnsafeHandler):
    def get(self, *args, **kwargs):
        current_version = self.dataset.current_version()

        ret = {
            'short_name': self.dataset.short_name,
            'full_name': self.dataset.full_name,
            'description': current_version.description,
            'terms': current_version.terms
        }

        self.finish(json.dumps(ret))


class GetUser(handlers.UnsafeHandler):
    def get(self, *args, **kwargs):
        user = self.current_user

        ret = {
                'user': None,
                'email': None,
                'trusted': False,
                'admin': False,
                'has_requested_access': False
        }
        if user:
            ### TODO there should probably be another way to figure out whether
            ## someone already has access or not. REST-endpoint or something
            ## similar, not really sure yet how this should be handled. I'm adding
            ## it here now so we can get the information to the browser.

            has_requested_access = False
            try:
                db.DatasetAccess.select().where(
                        db.DatasetAccess.user == user,
                        db.DatasetAccess.dataset == self.dataset).get()
                has_requested_access = True
            except:
                has_requested_access = False


            ret = {
                    'user':         user.name,
                    'email':        user.email,
                    'trusted':      self.is_authorized(),
                    'admin':        self.is_admin(),
                    'has_requested_access': has_requested_access
            }

        logging.info("getUser: " + str(ret['user']) + ' ' + str(ret['email']))
        self.finish(json.dumps(ret))

class CountryList(handlers.UnsafeHandler):
    def get(self, *args, **kwargs):
        self.write({'countries': [{'name': c} for c in self.country_list()]})

    def country_list(self):
        return ["Afghanistan", "Albania", "Algeria", "American Samoa", "Andorra",
                "Angola", "Anguilla", "Antarctica", "Antigua and Barbuda",
                "Argentina", "Armenia", "Aruba", "Australia", "Austria",
                "Azerbaijan", "Bahamas", "Bahrain", "Bangladesh", "Barbados",
                "Belarus", "Belgium", "Belize", "Benin", "Bermuda", "Bhutan",
                "Bolivia", "Bosnia and Herzegovina", "Botswana", "Brazil",
                "British Indian Ocean Territory", "British Virgin Islands",
                "Brunei", "Bulgaria", "Burkina Faso", "Burundi", "Cambodia",
                "Cameroon", "Canada", "Cape Verde", "Cayman Islands",
                "Central African Republic", "Chad", "Chile", "China",
                "Christmas Island", "Cocos Islands", "Colombia", "Comoros",
                "Cook Islands", "Costa Rica", "Croatia", "Cuba", "Curacao",
                "Cyprus", "Czech Republic", "Democratic Republic of the Congo",
                "Denmark", "Djibouti", "Dominica", "Dominican Republic",
                "East Timor", "Ecuador", "Egypt", "El Salvador",
                "Equatorial Guinea", "Eritrea", "Estonia", "Ethiopia",
                "Falkland Islands", "Faroe Islands", "Fiji", "Finland", "France",
                "French Polynesia", "Gabon", "Gambia", "Georgia", "Germany",
                "Ghana", "Gibraltar", "Greece", "Greenland", "Grenada", "Guam",
                "Guatemala", "Guernsey", "Guinea", "Guinea-Bissau", "Guyana",
                "Haiti", "Honduras", "Hong Kong", "Hungary", "Iceland", "India",
                "Indonesia", "Iran", "Iraq", "Ireland", "Isle of Man", "Israel",
                "Italy", "Ivory Coast", "Jamaica", "Japan", "Jersey", "Jordan",
                "Kazakhstan", "Kenya", "Kiribati", "Kosovo", "Kuwait",
                "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", "Liberia",
                "Libya", "Liechtenstein", "Lithuania", "Luxembourg", "Macau",
                "Macedonia", "Madagascar", "Malawi", "Malaysia", "Maldives",
                "Mali", "Malta", "Marshall Islands", "Mauritania", "Mauritius",
                "Mayotte", "Mexico", "Micronesia", "Moldova", "Monaco",
                "Mongolia", "Montenegro", "Montserrat", "Morocco", "Mozambique",
                "Myanmar", "Namibia", "Nauru", "Nepal", "Netherlands",
                "Netherlands Antilles", "New Caledonia", "New Zealand",
                "Nicaragua", "Niger", "Nigeria", "Niue", "North Korea",
                "Northern Mariana Islands", "Norway", "Oman", "Pakistan", "Palau",
                "Palestine", "Panama", "Papua New Guinea", "Paraguay", "Peru",
                "Philippines", "Pitcairn", "Poland", "Portugal", "Puerto Rico",
                "Qatar", "Republic of the Congo", "Reunion", "Romania", "Russia",
                "Rwanda", "Saint Barthelemy", "Saint Helena",
                "Saint Kitts and Nevis", "Saint Lucia", "Saint Martin",
                "Saint Pierre and Miquelon",
                "Saint Vincent and the Grenadines", "Samoa", "San Marino",
                "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia",
                "Seychelles", "Sierra Leone", "Singapore", "Sint Maarten",
                "Slovakia", "Slovenia", "Solomon Islands", "Somalia",
                "South Africa", "South Korea", "South Sudan", "Spain",
                "Sri Lanka", "Sudan", "Suriname", "Svalbard and Jan Mayen",
                "Swaziland", "Sweden", "Switzerland", "Syria", "Taiwan",
                "Tajikistan", "Tanzania", "Thailand", "Togo", "Tokelau", "Tonga",
                "Trinidad and Tobago", "Tunisia", "Turkey", "Turkmenistan",
                "Turks and Caicos Islands", "Tuvalu", "U.S. Virgin Islands",
                "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom",
                "United States", "Uruguay", "Uzbekistan", "Vanuatu", "Vatican",
                "Venezuela", "Vietnam", "Wallis and Futuna", "Western Sahara",
                "Yemen", "Zambia", "Zimbabwe" ];


class RequestAccess(handlers.SafeHandler):
    def get(self, *args, **kwargs):
        user = self.current_user
        name = user.name
        email = user.email

        logging.info("Request: " + name + ' ' + email)
        self.finish(json.dumps({'user':name, 'email':email}))

    def post(self, *args, **kwargs):
        userName    = self.get_argument("userName", default='',strip=False)
        email       = self.get_argument("email", default='', strip=False)
        affiliation = self.get_argument("affiliation", strip=False)
        country     = self.get_argument("country", strip=False)
        newsletter  = self.get_argument("newsletter", strip=False)

        # This is the only chance for XSRF in the application
        # avoid it by checking that the email sent by google is the same as
        # supplied by the form post
        user = self.current_user
        if user.email != email:
            return

        user.affiliation = affiliation
        user.country = country
        logging.info("Inserting into database: {}, {}".format(user.name, user.email))

        try:
            with db.database.atomic():
                user.save() # Save to database
                db.DatasetAccess.create(
                        user             = user,
                        dataset          = self.dataset,
                        wants_newsletter = newsletter
                    )
                db.UserLog.create(
                        user = user,
                        dataset = self.dataset,
                        action = 'access_requested'
                    )
        except Exception as e:
            logging.error(e)


class LogEvent(handlers.SafeHandler):
    def get(self, sEvent):
        user = self.current_user

        ok_events = ['download','consent']
        if sEvent in ok_events:
            db.UserLog.create(
                    user = user,
                    dataset = self.dataset,
                    action = sEvent
                )
        else:
            raise tornado.web.HTTPError(400, reason="Can't log that")

class ApproveUser(handlers.AdminHandler):
    def get(self, sEmail):
        with db.database.atomic():
            user = db.User.select().where(db.User.email == sEmail).get()

            da = db.DatasetAccess.select().where(
                        db.DatasetAccess.user == user,
                        db.DatasetAccess.dataset == self.dataset
                ).get()
            da.has_access = True
            da.save()

            db.UserLog.create(
                    user = user,
                    dataset = self.dataset,
                    action = 'access_granted'
                )

        msg = email.mime.multipart.MIMEMultipart()
        msg['to'] = sEmail
        msg['from'] = settings.from_address
        msg['subject'] = 'Swefreq account created'
        msg.add_header('reply-to', settings.reply_to_address)
        body = "Your Swefreq account has been activated."
        msg.attach(MIMEText(body, 'plain'))

        server = smtplib.SMTP(settings.mail_server)
        server.sendmail(msg['from'], [msg['to']], msg.as_string())


class RevokeUser(handlers.AdminHandler):
    def get(self, sEmail):
        if self.current_user.email == sEmail:
            # Don't let the admin delete hens own account
            return

        with db.database.atomic():
            user = db.User.select().where(db.User.email == sEmail).get()

            da = db.DatasetAccess.select(
                    ).join(
                        db.User
                    ).where(
                        db.User.email == sEmail,
                        db.DatasetAccess.dataset == self.dataset
                    ).get()
            da.delete_instance()

            db.UserLog.create(
                    user = user,
                    dataset = self.dataset,
                    action = 'access_revoked'
                )

class GetOutstandingRequests(handlers.SafeHandler):
    def get(self, *args, **kwargs):
        requests = db.get_outstanding_requests(self.dataset)

        json_response = []
        for request in requests:
            apply_date = request.apply_date.strftime('%Y-%m-%d')
            json_response.append({
                'user':        request.name,
                'email':       request.email,
                'affiliation': request.affiliation,
                'country':     request.country,
                'applyDate':   apply_date
            })

        self.finish(json.dumps(json_response))

class GetApprovedUsers(handlers.SafeHandler):
    def get(self, *args, **kwargs):
        ## All users that have access to the dataset and how many times they have
        ## downloaded it
        query = db.User.select(
                db.User, db.DatasetAccess.wants_newsletter
            ).join(
                db.DatasetAccess
            ).switch(
                db.User
            ).join(
                db.UserLog,
                peewee.JOIN.LEFT_OUTER,
                on=(   (db.User.user        == db.UserLog.user)
                     & (db.UserLog.action   == 'download')
                     & (db.UserLog.dataset  == db.DatasetAccess.dataset)
                )
            ).where(
                db.DatasetAccess.dataset    == self.dataset,
                db.DatasetAccess.has_access == 1
            ).annotate(db.UserLog)

        json_response = []
        for user in query:
            json_response.append({
                    'user':          user.name,
                    'email':         user.email,
                    'affiliation':   user.affiliation,
                    'country':       user.country,
                    'downloadCount': user.count,
                    'newsletter':    user.dataset_access.wants_newsletter
                })

        self.finish(json.dumps(json_response))
