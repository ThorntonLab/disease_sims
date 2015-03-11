#ifndef __TFL_FITNESS_MODELS_HPP__
#define __TFL_FITNESS_MODELS_HPP__

enum class MODEL { GENE_RECESSIVE = 0, GENE_ADDITIVE, MULTIPLICATIVE, POPGEN, EYREWALKER };

#include <diseaseSims/mutation_with_age.hpp>
#include <gene_based_model.hpp>
#include <multiplicative_model.hpp>
#include <popgen_model.hpp>
#include <EyreWalker_fitness.hpp>
#endif
