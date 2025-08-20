#ifndef GUARD_gen_models_h
#define GUARD_gen_models_h

#include "structs.h"

class GenModel {
    public:
        // note: lambda is the linear predictor, theta is the model parameter
        // i.e.: y ~ dist(theta), theta = inv_link(lambda)
        virtual double link(double param) = 0;
        virtual double inv_link(double lambda) = 0;
        virtual double log_lik(double y, double lambda) = 0;
        virtual double score(double y, double theta) = 0;
        virtual double jacobian(double theta) = 0;
  
        double proposal_mu_single(double &m, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi);
        double proposal_var_single(double m, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi);
        double proposal_mu_multi(double &m, const int &nid, suff_stat &ss, int &ensm_id, data_info &di, tree_prior_info &tree_pi);
        double proposal_var_multi(double m, const int &nid, suff_stat &ss, int &ensm_id, data_info &di, tree_prior_info &tree_pi);
        // double compute_log_lik(double &mu, suff_stat &ss, int &nid, data_info &di, tree_prior_info &tree_pi);
        double compute_node_lik_single(double &mu, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi);
    };

class Logit : public GenModel {
    public:
        double score(double y, double theta);
        double jacobian(double theta);
        double link(double param);
        double inv_link(double lambda);
        double log_lik(double y, double lambda);
    };

#endif