(function() {
    angular.module("App")
    .factory("DatasetVersions", function($http) {
        return {
            getDatasetVersions: getDatasetVersions,
        };
        function getDatasetVersions(dataset) {
            return $http.get("/api/datasets/" + dataset + "/versions")
                .then(function(data) {
                    return data.data.data;
            });
        };
    });
})();
