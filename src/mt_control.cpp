#include <ea/digital_evolution.h>
#include <ea/cmdline_interface.h>
#include <ea/subpopulation_founder.h>
#include <ea/line_of_descent.h>
#include <ea/analysis/archive.h>
#include <ea/metapopulation.h>


#include "multi_birth_selfrep_not_remote_ancestor.h"
#include "mt.h"
#include "gls.h"

using namespace ealib;


//! Configuration object for an EA.
struct lifecycle : public default_lifecycle {

    //! Called as the final step of EA construction (must not depend on configuration parameters)
    template <typename EA>
    void after_initialization(EA& ea) {
        if(ea.isa().size()) {
            return;
        }

        /* These are the instructions available. The ones common to most of Avida are defined in:
         ealib/libea/include/ea/digital_evolution/instruction_set.h */
        using namespace instructions;
        append_isa<nop_a>(0,ea);
        append_isa<nop_b>(0,ea);
        append_isa<nop_c>(0,ea);
        append_isa<nop_x>(ea);
        append_isa<mov_head>(ea);
        append_isa<if_label>(ea);
        append_isa<h_search>(ea);
        append_isa<nand>(ea);
        append_isa<push>(ea);
        append_isa<pop>(ea);
        append_isa<swap>(ea);
        append_isa<inc>(ea);
        append_isa<dec>(ea);
        // communication
        // send message
        append_isa<tx_msg>(ea);
        append_isa<tx_msg_check_task>(ea);
        // receive message
        append_isa<rx_msg>(ea);
        // broadcast message
        append_isa<bc_msg>(ea);
        append_isa<bc_msg_check_task>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea);
        append_isa<h_alloc>(ea);
        append_isa<h_copy>(ea);
        // divide producing a cell
        append_isa<h_divide_local>(ea);
        // take in inputs from the environment -- fixed inputs means the same numbers.
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);

        // used for multicell / organism replication
        append_isa<donate_res_to_group>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<is_neighbor>(ea);
        /* Divide producing another organism. For this to happen the organism (or group) must
         have acquired XXX resources. How this works: the cells perform tasks. Tasks gain resources.
         These reosources belong to the cell until they donate them to the organism (group) see
         above instruction donate_res_to_group. */
        append_isa<h_divide_remote>(ea);
        append_isa<become_soma>(ea);
        append_isa<if_germ>(ea);
        append_isa<if_soma>(ea);
        append_isa<if_res_more_than_thresh>(ea);
        append_isa<if_res_less_than_thresh>(ea);

        add_event<task_resource_consumption>(ea);
        add_event<task_switching_cost>(ea);


        add_event<ts_birth_event>(ea);
        add_event<task_mutagenesis>(ea); // <- this event mutates genomes when tasks are performed
        add_event<gs_inherit_event>(ea); // <- this event controls the germ/soma state being inherited

        typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
        typedef typename EA::resource_ptr_type resource_ptr_type;

        // Add tasks -- all of the additive values are set to 0 since individual cell's don't change their
        // execution speed. We are just using tasks to give the cells resources that can be used
        // for replication.
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<0> >("and", ea);
        task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<0> >("ornot", ea);
        task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<0> >("or", ea);
        task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<0> >("andnot", ea);
        task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<0> >("nor", ea);
        task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<0> >("xor", ea);
        task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<0> >("equals", ea);

        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
        double init_amt = get<RES_INITIAL_AMOUNT>(ea, 0);
        double inflow = get<RES_INFLOW_AMOUNT>(ea,0);
        double outflow = get<RES_OUTFLOW_FRACTION>(ea,0);
        double frac = get<RES_FRACTION_CONSUMED>(ea,0);

        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
        resource_ptr_type resA = make_resource("resA", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resB = make_resource("resB", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resC = make_resource("resC", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resD = make_resource("resD", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resE = make_resource("resE", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resF = make_resource("resF", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resG = make_resource("resG", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resH = make_resource("resH", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resI = make_resource("resI", init_amt, inflow, outflow, frac, ea);

        // Set the mutagenic effects of each task.
        put<TASK_MUTATION_MULT>(get<NOT_MUTATION_MULT>(ea), *task_not);
        put<TASK_MUTATION_MULT>(get<NAND_MUTATION_MULT>(ea), *task_nand);
        put<TASK_MUTATION_MULT>(get<AND_MUTATION_MULT>(ea), *task_and);
        put<TASK_MUTATION_MULT>(get<ORNOT_MUTATION_MULT>(ea), *task_ornot);
        put<TASK_MUTATION_MULT>(get<OR_MUTATION_MULT>(ea), *task_or);
        put<TASK_MUTATION_MULT>(get<ANDNOT_MUTATION_MULT>(ea), *task_andnot);
        put<TASK_MUTATION_MULT>(get<NOR_MUTATION_MULT>(ea), *task_nor);
        put<TASK_MUTATION_MULT>(get<XOR_MUTATION_MULT>(ea), *task_xor);
        put<TASK_MUTATION_MULT>(get<EQUALS_MUTATION_MULT>(ea), *task_equals);

        task_not->consumes(resA);
        task_nand->consumes(resB);
        task_and->consumes(resC);
        task_ornot->consumes(resD);
        task_or->consumes(resE);
        task_andnot->consumes(resF);
        task_nor->consumes(resG);
        task_xor->consumes(resH);
        task_equals->consumes(resI);

    }

};


/* This code is necessary for line of descent tracking. */
template <typename T>
struct subpop_trait : subpopulation_founder_trait<T>, fitness_trait<T> {
    typedef subpopulation_founder_trait<T> parent1_type;
    typedef fitness_trait<T> parent2_type;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp("subpopulation_founder_trait", boost::serialization::base_object<parent1_type>(*this));
        ar & boost::serialization::make_nvp("fitness_trait", boost::serialization::base_object<parent2_type>(*this));
    }
};



/* This defines the organism, which is a group of cells. */
typedef digital_evolution
< lifecycle
, recombination::asexual
, round_robin
, multibirth_selfrep_not_remote_ancestor // the ancestor does one task repeatedly until it replicate and form another organism. You can find it in multibirth_selfrep_not_remote_ancestor.h

// Kate - this controls whether a cell can replicate over another cell
, empty_facing_neighbor   // use for control
//, faced_neighbor   // use for rogue cell treatment

, generate_single_ancestor
> sea_type;

/* This defines our population of organisms. */
typedef metapopulation
< sea_type
, quiet_nan
, mutation::operators::no_mutation
, quiet_nan
, generational_models::isolated_subpopulations
, ancestors::default_subpopulation
, dont_stop
, fill_metapopulation
, default_lifecycle
, subpop_trait
> mea_type;




/*!
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        add_option<SPATIAL_X>(this); // number of cells an organism is in x dimension (usually 5)
        add_option<SPATIAL_Y>(this); // number of cells an organism is in y dimension (usually 5)
        add_option<METAPOPULATION_SIZE>(this); // number or organisms (usually 1000)
        add_option<POPULATION_SIZE>(this); // max number of cells in an organism (usually 25)
        add_option<REPRESENTATION_SIZE>(this);
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<SCHEDULER_RESOURCE_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this); // used to configure probability of mutation during cell replication
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<MUTATION_UNIFORM_INT_MIN>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this); // must be the number of instructions used. Otherwise if > errors.

        add_option<ANALYSIS_INPUT>(this);


        // ts specific options
        add_option<TASK_SWITCHING_COST>(this); // set to 0 for mt projects, but useful if you want to study task switching
        add_option<GERM_MUTATION_PER_SITE_P>(this); // the mutation rate used for organismal replication
        add_option<GROUP_REP_THRESHOLD>(this); // the number of resources required for a group (really an organism) to replicate

        // Austin - this is the config that when raised high enough results in multicells reverting to unicells
        add_option<IND_REP_THRESHOLD>(this); // the number of resources required for an individual (cell) to replicate

        // dirty work mutation stuff
        add_option<TASK_MUTATION_PER_SITE_P>(this); // the base rate of mutations for a task
        // each task has it's own mutation rate which is the base rate * its multiplication rate
        // so for example, with base rate b, task NOT's rate would be NOT_MUTATION_MULT * b
        add_option<NOT_MUTATION_MULT>(this); // normally set to 0
        add_option<NAND_MUTATION_MULT>(this);
        add_option<AND_MUTATION_MULT>(this);
        add_option<ORNOT_MUTATION_MULT>(this);
        add_option<OR_MUTATION_MULT>(this);
        add_option<ANDNOT_MUTATION_MULT>(this);
        add_option<NOR_MUTATION_MULT>(this);
        add_option<XOR_MUTATION_MULT>(this);
        add_option<EQUALS_MUTATION_MULT>(this);

        add_option<RES_INITIAL_AMOUNT>(this);
        add_option<RES_INFLOW_AMOUNT>(this);
        add_option<RES_OUTFLOW_FRACTION>(this);
        add_option<RES_FRACTION_CONSUMED>(this);
        add_option<COST_RAMP>(this);
        add_option<COST_START_UPDATE>(this);

        add_option<ARCHIVE_INPUT>(this);
        add_option<ARCHIVE_OUTPUT>(this);
    }

    virtual void gather_tools() {
        // these are examples for when you need tools - they can be used to reload populations, reload
        // line of descents and then perform further analyses.

        //add_tool<movie>(this);
        //add_tool<ealib::analysis::lod_knockouts>(this);


    }

    virtual void gather_events(EA& ea) {
        add_event<mt_gls_propagule>(ea); // this event defines how organism replication occurs.

        // this is used for line of descent tracking. It's memory expensive. Normally I run with 1 gig. If
        // I use LOD - it ups it to 4 gigs.
//        add_event<datafiles::mrca_lineage>(ea);
//        add_event<subpopulation_founder_event>(ea);
        add_event<task_performed_tracking>(ea);

    }
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
