(function() {
    angular.module("App")
    .controller("browserController", ["$routeParams", "User", "Dataset", "browserTranscript", "browserRegion",
                             function( $routeParams,   User,   Dataset,   browserTranscript,   browserRegion) {
        var localThis = this;
        localThis.browserLink = browserLink;
        localThis.dataset = {'shortName':$routeParams.dataset};
        activate();

        function activate() {
            User.getUser().then(function(data) {
                localThis.user = data;
            });
            if ($routeParams.transcript) {
                browserTranscript.getTranscript($routeParams.dataset, $routeParams.transcript).then( function(data) {
                    localThis.transcript = data.transcript;
                    localThis.gene       = data.gene;
                });
            }
            if ($routeParams.region) {
                browserRegion.getRegion($routeParams.dataset, $routeParams.region).then( function(data) {
                    localThis.region = data.region;
                });
            }
            Dataset.getDataset($routeParams.dataset, $routeParams.version)
                .then(function(data) {
                    localThis.dataset = data.dataset;
                },
                function(error) {
                    localThis.error = error;
            });
        }

        function browserLink(link) {
            return "/dataset/" + $routeParams.dataset + "/browser/" + link;
        }
    }]);
})();
