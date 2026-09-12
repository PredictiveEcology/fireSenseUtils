## The National Fire Database point archive (all fires, shapefile) on the CFS server.
##
## The file was `NFDB_point.zip`; CFS renamed it `NFDB_point_shp.zip` (with `_txt` and `large_fires`
## variants beside it), and the old name returns HTTP 404. The re-download check in
## getFirePoints_NFDB() and getFirePoints_NFDB_V2() then could not fetch a newer release, so the
## last copy on disk (fires to 2024) kept being used.
nfdbPointUrl <- function() {
  "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point_shp.zip"
}
