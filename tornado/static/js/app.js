(function() {

    /////////////////////////////////////////////////////////////////////////////////////
    // create the module and name it App
    var App = angular.module('App', ['ngRoute', 'ngCookies'])
    var gData = {'userName':'',
                 'email':'',
                 'affiliation':'',
                 'trusted':false,
                 'admin':false}

    /////////////////////////////////////////////////////////////////////////////////////
    App.directive('consent', function ($cookies) {
        return {
            scope: {},
            template:
                '<div style="position: relative; z-index: 1000">' +
                '<div style="background: #eee; position: fixed; bottom: 0; left: 0; right: 0; height: 20px" ng-hide="consent()">' +
                '<span style="margin-left: 5px;">This site uses cookies, please see our <a href="https://swefreq.nbis.se/#/privacyPolicy/">privacy policy</a>, ok, I <a href="" ng-click="consent(true)">agree.</a></span>' +
                '<span style="float: right;margin-right: 5px;"><a href="" ng-click="consent(true)">X</a></span>' +
                '</div>' +
                '</div>',
            controller: function ($scope) {
                var _consent = $cookies.get('consent');
                $scope.consent = function (consent) {
                    if (consent === undefined) {
                        return _consent;
                    } else if (consent) {
                        $cookies.put('consent', true);
                        _consent = true;
                    }
                };
            }
        };
    });

    // Modified from
    // http://stackoverflow.com/questions/16199418/how-to-set-bootstrap-navbar-active-class-with-angular-js
    App.directive('bsActiveLink', ['$location', function($location) {
        return {
            restrict: 'A',
            replace: false,
            link: function (scope, elem) {
                // After the route has changed
                scope.$on("$routeChangeSuccess", function () {
                    $location.path()
                    var hrefs = ['/#' + $location.path(),
                                 '#' + $location.path(), //html5: false
                                 $location.path()]; //html5: true

                    angular.forEach(elem.find('a'), function(a) {
                        a = angular.element(a);
                        console.log(a.attr('href'));
                        if ( -1 != hrefs.indexOf(a.attr('href')) ) {
                            a.parent().addClass('active');
                        } else {
                            a.parent().removeClass('active');
                        }
                    });
                })
            }
        }
    }]);


    App.controller('mainController', function($http, $scope) {
        var localThis = this;
        localThis.data = gData;

        this.getUsers = function(){
            $http.get('/api/users/me').success(function(data){
                console.log(data);
                localThis.data.userName = data.user;
                localThis.data.email = data.email;
                localThis.data.trusted = data.trusted;
                localThis.data.has_requested_access = data.has_requested_access;
                localThis.data.admin = data.admin;
            });
        };
        this.getUsers();
    });

    /////////////////////////////////////////////////////////////////////////////////////

    App.controller('homeController', function($http, $scope, $sce) {
        var localThis = this;
        localThis.datasets = [];
        localThis.getDatasets = function(){
            $http.get('/api/datasets').success(function(res){
                var len = res.data.length;
                for (var i = 0; i < len; i++) {
                    d = res.data[i];
                    d.version.description = $sce.trustAsHtml(d.version.description)

                    localThis.datasets.push(d);
                }
            });
        };
        localThis.getDatasets();
    });

    /////////////////////////////////////////////////////////////////////////////////////

    App.controller('adminController', function($http, $scope) {
        var localThis = this;
        this.userName = '';
        this.email = '';
        localThis.data = gData;

        this.getUsers = function(){
            $http.get('/api/users/me').success(function(data){
                console.log(data);
                localThis.data.userName = data.user;
                localThis.data.email = data.email;
                localThis.data.trusted = data.trusted;
                localThis.data.has_requested_access = data.has_requested_access;
                localThis.data.admin = data.admin;
                if(data.admin == true){
                    // TODO: Change this to one call that is then filtered into
                    // the two different datasets? Or just filter in the view.
                    // This is currently broken.
                    $http.get('/api/datasets/swegen/users').success(function(data){
                        localThis.data.requests = data;
                    });
                    $http.get('/api/datasets/swegen/users').success(function(data){
                        localThis.data.approvedUsers = data;
                        localThis.data.emails = []
                        for (var idx in data) {
                            var user = data[idx];
                            if (user.newsletter == 1) {
                                localThis.data.emails.push(user['email']);
                            }
                        }
                    });
                };
            });
        };
        this.getUsers();

        this.revokeUser = function(userData){
            $http.get('/api/datasets/swegen/users/' + userData.email + '/revoke').success(function(data){
                localThis.getUsers();
            });
        };

        this.approvedUser = function(userData){
            $http.get('/api/datasets/swegen/users/' + userData.email + '/approve').success(function(data){
                $http.get('/api/datasets/swegen/users/').success(function(data){
                    localThis.getUsers();
                });
            });
        };
    });

     /////////////////////////////////////////////////////////////////////////////////////

    App.controller('dataBeaconController', function($http, $window) {
        var beacon = this;
        beacon.pattern = { 'chromosome': "\\d+" };
        beacon.beacon_info = {};
        $http.get('/api/info').success(function(data) {
            beacon.beacon_info = data;
            beacon.datasets = data['datasets'];
            beacon.dataset = data['datasets'][0]['id'];
            beacon.reference = data['datasets'][0]['reference'];
        });
        beacon.search = function() {
            beacon.color = 'black';
            beacon.response = "Searching...";
            $http.get('/api/query', { 'params': { 'chrom': beacon.chromosome, 'pos': beacon.position - 1, 'allele': beacon.allele, 'referenceAllele': beacon.referenceAllele, 'dataset': beacon.dataset, 'ref': beacon.reference}})
                .then(function (response){
                    if (response.data['response']['exists']) {
                        beacon.response = "Present";
                        beacon.color = 'green';
                    }
                    else {
                        beacon.response = "Absent";
                        beacon.color = "red";
                    }
                },
                function (response){
                    beacon.response="Error";
                    beacon.color = 'black';
                });
        }
    });

     /////////////////////////////////////////////////////////////////////////////////////

    App.controller('downloadDataController', function($http, $scope, $sce) {
        this.lChecked = true;
        var localThis = this;
        localThis.data = gData;
        this.isChecked = function(){
            if(localThis.lChecked){
                $http.get('/api/log/consent').success(function(data){
                    console.log('Consented');
                });
            }
            localThis.lChecked = false;
            localThis.checked = true;
        };

        this.downloadData = function(){
            $http.get('/api/log/download').success(function(data){
                console.log("Downloading")
            });
        };

        localThis.getDataset = function(){
            $http.get('/api/datasets/swegen').success(function(data){
                localThis.short_name  = data.short_name;
                localThis.full_name   = data.full_name;
                localThis.description = data.description;
                localThis.terms       = data.terms;
                localThis.version     = data.version;
                localThis.has_image   = data.has_image;
                localThis.files = data.files;
            });
        };
        localThis.getDataset();
        localThis.getTerms = function(){
            return $sce.trustAsHtml(localThis.terms);
        };
    });

    /////////////////////////////////////////////////////////////////////////////////////

    App.controller('requestController', function($http, $scope, $location) {
        var localThis = this;
        localThis.data = gData;
        localThis.data.newsletter = true;
        $http.get('/api/countries').success(function(data) {
            localThis.data['availableCountries'] = data['countries'];
        });

        this.sendRequest = function(valid){
            if (!valid) {
                return;
            }
            $http.defaults.headers.post["Content-Type"] = "application/x-www-form-urlencoded";
            $http({url:'/api/datasets/swegen/users/' + localThis.data.email + '/request',
                   method:'POST',
                   data:$.param({'email':localThis.data.email,
                                 'userName':localThis.data.userName,
                                 'affiliation':localThis.data.affiliation,
                                 'country': localThis.data.country['name'],
                                 'newsletter': localThis.data.newsletter ? 1 : 0
                        })
                })
                .success(function(data){
                    console.log(data);
                    $location.path("/addedRequest");
                });
        };
    });

    /////////////////////////////////////////////////////////////////////////////////////

    App.controller('addedRequestController', function($http, $scope) {
        var localThis = this;
    });

    ////////////////////////////////////////////////////////////////////////////
    // configure routes
    App.config(function($routeProvider, $locationProvider) {
        $routeProvider
            .when('/',               { templateUrl: 'static/js/ng-templates/home.html'          })
            .when('/dataBeacon/',    { templateUrl: 'static/js/ng-templates/dataBeacon.html'    })
            .when('/downloadData/',  { templateUrl: 'static/js/ng-templates/downloadData.html'  })
            .when('/requestAccess/', { templateUrl: 'static/js/ng-templates/requestAccess.html' })
            .when('/addedRequest/',  { templateUrl: 'static/js/ng-templates/addedRequest.html'  })
            .when('/privacyPolicy/', { templateUrl: 'static/js/ng-templates/privacyPolicy.html' })
            .when('/admin/',         { templateUrl: 'static/js/ng-templates/admin.html'         })
            .when('/terms/',         { templateUrl: 'static/js/ng-templates/terms.html'         })
            .otherwise(              { templateUrl: 'static/js/ng-templates/404.html'           });

        // Use the HTML5 History API
        $locationProvider.html5Mode(true);
    });
})();

