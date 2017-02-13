index="""
<!DOCTYPE html>
<!-- define angular app -->
<html ng-app="App">
  <head>
    <title>SweFreq</title>
    <!-- SCROLLS -->
    <meta charset="utf-8" />
    <!-- Bootstrap -->
    <link rel="stylesheet" href="//netdna.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css" />
    <link rel="stylesheet" href="//netdna.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.css" />
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js"></script>
    <!-- Angular -->
    <script src="https://ajax.googleapis.com/ajax/libs/angularjs/1.5.0/angular.min.js"></script>
    <script src="https://ajax.googleapis.com/ajax/libs/angularjs/1.5.0/angular-route.min.js"></script>
    <script src="https://ajax.googleapis.com/ajax/libs/angularjs/1.5.0/angular-cookies.min.js"></script>
    <!-- Angular grid gui -->
    <link rel="stylesheet" href="/javascript/main.css" type="text/css" />
    <link rel="stylesheet" href="/javascript/local.css" type="text/css" />
    <!-- The application -->
    <script src="javascript/app.js"></script>

    <!-- Google Analytics -->
    <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
        (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
        m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
    ga('create', 'UA-85976442-1', 'auto');
    ga('send', 'pageview');
    </script>
    <!-- End Google Analytics -->
  </head>
  <!-- define angular controller -->
  <body ng-controller="mainController as mainCtrl">
    <nav class="navbar navbar-default">
      <div class="container">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed"
          data-toggle="collapse" data-target="#navbar-things"
          aria-expanded="false">
              <span class="sr-only">Toggle navigation</span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="/#/"><img height="50" src="/javascript/swefreq_ed_50px.png"></a>
        </div>

        <div class="collapse navbar-collapse" id="navbar-things">
        <ul class="nav navbar-nav navbar-right" bs-active-link>
          <li><a href="/#/about/">About</a>
          <li><a href="/#/terms/">Terms of use</a>
          <li><a href="/#/dataBeacon/">Data Beacon</a>
          <li><a href="{{ExAC}}">SweGen Browser</a>
        {% if has_access %}
          <li><a href="/#/downloadData/">Download Data</a>
        {% elif user_name %}
          <li><a href="/#/requestAccess/">Request Access</a>
        {% end %}
        {% if is_admin %}
          <li><a href="/#/admin/">Admin</a>
        {% end %}
        {% if user_name %}
          <li><a href="/logout" title="Logout {{user_name}}">Logout</a>
        {% else %}
          <li><a href="/login">Login</a>
        {% end %}
        </ul>
        </div>
      </div>
    </nav>
    <div id="main">
      <!-- angular templating -->
      <!-- this is where content will be injected by angular -->
      <div ng-view=""></div>
    </div>
      </div>
    <footer class="text-center"></footer>
  </body>
</html>
"""

