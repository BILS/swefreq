from . import browser_handlers as handlers

# Browser links
routes = [(r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/gene/(?P<gene>[^/]+)" ,                              handlers.GetGene),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/region/(?P<region>[^\/]+)",                          handlers.GetRegion),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/transcript/(?P<transcript>[^/]+)",                   handlers.GetTranscript),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/variant/(?P<variant>[^/]+)",                         handlers.GetVariant),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/variants/(?P<datatype>[^/]+)/(?P<item>[^/]+)",       handlers.GetVariants),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/coverage/(?P<datatype>[^/]+)/(?P<item>[^/]+)",       handlers.GetCoverage),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/coverage_pos/(?P<datatype>[^/]+)/(?P<item>[^/]+)",   handlers.GetCoveragePos),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/autocomplete/(?P<query>[^/]+)",                      handlers.Autocomplete),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/search/(?P<query>[^/]+)",                            handlers.Search),
          (r"/api/datasets/(?P<dataset>[^/]+)/(?:version/(?P<ds_version>[^/]+)/)?browser/download/(?P<datatype>[^/]+)/(?P<item>[^/]+)",       handlers.Download),
        ]
