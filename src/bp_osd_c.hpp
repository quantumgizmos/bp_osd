#ifndef BP_OSD_C_H
#define BP_OSD_C_H

#include <json.hpp>
#include "timing.h"
#include "json_util.h"
using json = nlohmann::json;

//C include
extern "C" {
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "load_alist.h"
#include "syndrome.h"
#include "binary_char.h"
#include "bp_decoder_ms.h"
#include "mtwister.h"
#include "mod2sparse_extra.h"
#include "logical_check.h"
}

#include "bp_osd.h"
#include "setup_mtwister.h"
#include "sim_util.h"

#endif //BP_OSD_C_H
