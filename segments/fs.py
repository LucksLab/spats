INVALID = 0
error = 13
Indel = 1383
mutation = 1384
insertion = 1385
deletion = 1386
substitution = 1387
NGSMutation = 1388
NGSRead = 1389
goner = 1390

g_fsMap = {
    goner : "goner",
    error : "error",
    Indel : "Indel",
    mutation : "mutation",
    insertion : "insertion",
    deletion : "deletion",
    substitution : "substitution",
    NGSMutation : "NGSMutation",
    NGSRead : "NGSRead",
    goner : "goner",
}

def toString(fs):
    return g_fsMap.get(fs, "<INVALID>")

def fromString(arg):
    return globals()[arg]
