(function() {
    angular.module("App")
    .factory("DatasetFiles", ["$http", function($http) {
        return {
            getFiles: getFiles,
        };

        function getFiles(dataset, version) {
            var fileUri = "/api/dataset/" + dataset + "/files";
            if ( version ) {
                fileUri = "/api/dataset/" + dataset + "/versions/" + version + "/files";
            }
            return $http.get(fileUri).then(function(data) {
                return data.data.files;
            });
        }
    }]);
})();
