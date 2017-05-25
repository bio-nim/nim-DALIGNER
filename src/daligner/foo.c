#include "DB.h"
#include <string.h>

int main(){
    HITS_DB db;
    db.loaded = db.bases = db.reads = db.tracks = 0;
    char const* name0 = "tmpdb.db";
    char* name = strdup(name0);
    Open_DB(name, &db);
    Close_DB(&db);
    return 0;
}
