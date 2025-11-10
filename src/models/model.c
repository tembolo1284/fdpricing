/**
 * model.c - Common model functions
 */

#include "fdpricing.h"
#include "internal/models/model.h"
#include "internal/core/context.h"

void fdp_model_free(fdp_model_t* model)
{
    if (!model) return;
    
    /* Call model-specific destructor */
    if (model->vtable && model->vtable->destroy) {
        model->vtable->destroy(model);
    }
    
    /* Free the model structure itself */
    fdp_ctx_free(model->ctx, model);
}

fdp_model_type_t fdp_model_get_type(const fdp_model_t* model)
{
    return model ? model->type : FDP_MODEL_GBM;
}

int fdp_model_is_1d(const fdp_model_t* model)
{
    if (!model) return 0;
    
    return (model->type == FDP_MODEL_GBM ||
            model->type == FDP_MODEL_MERTON_JUMP ||
            model->type == FDP_MODEL_KOU_JUMP ||
            model->type == FDP_MODEL_LOCAL_VOL);
}

int fdp_model_is_2d(const fdp_model_t* model)
{
    if (!model) return 0;
    
    return (model->type == FDP_MODEL_HESTON ||
            model->type == FDP_MODEL_SABR ||
            model->type == FDP_MODEL_BATES);
}

int fdp_model_has_jumps(const fdp_model_t* model)
{
    if (!model) return 0;
    
    return (model->type == FDP_MODEL_MERTON_JUMP ||
            model->type == FDP_MODEL_KOU_JUMP ||
            model->type == FDP_MODEL_BATES);
}
