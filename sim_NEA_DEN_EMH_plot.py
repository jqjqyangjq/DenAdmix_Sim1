from collections import defaultdict, namedtuple
import msprime
import pandas as pd
import math
#import stdpopsim

#helper functions    from Leo
MERGE = namedtuple("merge", ('source', 'target', 'time'))

def define_popconfig(pop_params):
    """Generate list of population configurations for msprime."""
    demographic_events = msprime.Demography()
    for p in pop_params:
        demographic_events.add_population(name=str(p),initial_size=pop_params[p]["Ne"], initially_active=pop_params[p]["Init"])
    return demographic_events

def define_samples(samples):
    """Generate list of sample definitions for msprime."""
    sample_names = []
    for i, pop in enumerate(samples):
        times_generation = [t*1000//29 for t in samples[pop]["t_sample"]]
        sample_names.extend([msprime.SampleSet(samples[pop]["ind"],population=pop, time=int(t),ploidy=2) for t in times_generation])
    return sample_names


def demo_archaic_introgression(
    N_p,
    D_p,
    D_s,
    ):

    #population sizes, size changes and labels
    n_ea = 5054
    n_out_of_afr = 781
    n_afr = 27122
    n_archaic_d = 1500
    n_archaic_n = 1000
    pop_params = {
        "early_NEA":  {"id": 0, "Ne": n_archaic_n, "Init": True},
        "late_NEA": {"id": 1, "Ne": n_archaic_n, "Init": True },
        "Intro_NEA" : {"id": 2, "Ne": n_archaic_n, "Init": True},
        "late_DEN" : {"id": 3,"Ne": n_archaic_d, "Init": True},
        "early_DEN": {"id": 4, "Ne": n_archaic_d, "Init": True},
        "EA": {"id": 5, "Ne": n_ea, "Init": True},
        "AFR": {"id": 6, "Ne": n_afr, "Init": True },
        "Intro_DEN" : {"id": 7, "Ne": n_archaic_d, "Init": True},
        "Archaic": {"id": 8, "Ne": 10000, "Init": False},
        "Ancestral": {"id": 9, "Ne": 23275, "Init": False}
    }
    #ids = dict((k, p["id"]) for k, p in pop_params.items())

    demographic_events = define_popconfig(pop_params)

    """3. set up base tree"""
    base_tree = [
            MERGE('early_DEN','late_DEN', 250),
            MERGE('Intro_DEN', 'late_DEN', D_s),
            MERGE('early_NEA' ,'late_NEA', 150),
            MERGE('EA','AFR', 72),
            MERGE('Intro_NEA','late_NEA', 100)
            ]

    for m in base_tree:
        demographic_events.add_population_split(time=(m.time *1000//29), derived=[str(m.source)], ancestral=str(m.target))
    demographic_events.add_population_split(time=440000//29, derived=['late_NEA', 'late_DEN'], ancestral='Archaic')
    demographic_events.add_population_split(time=656908//29, derived=['AFR', 'Archaic'], ancestral='Ancestral')
    
    demographic_events.add_population_parameters_change(time=72000//29+100,initial_size=23275,population="AFR")
    #
    demographic_events.add_population_parameters_change(time=656908//29,initial_size=18296,population="Ancestral")
    
    """Out of Africa bottleneck"""
    demographic_events.add_population_parameters_change(time=72000//29-100,initial_size=781,population="EA")
    ###merge back wit afr, from base tree
    #demographic_events.add_population_parameters_change(time=72000//29,initial_size=23275,population="afr")
    """eruasian bottle neck"""
    demographic_events.add_population_parameters_change(time=57100//29,initial_size=7264,population="EA")
    
    demographic_events.add_population_parameters_change(time=57100//29-100,initial_size=1305,population="EA")
    
    demographic_events.add_mass_migration(time=(50000//29), source="EA", dest="Intro_NEA",proportion=N_p)
    demographic_events.add_mass_migration(time=(45000//29), source="EA", dest="Intro_DEN",proportion=D_p )



    demographic_events.sort_events()

    print("Printing demography")
    print(demographic_events.debug())
    check = demographic_events.debug()
    check.print_history(output=open('sim_demo_%s.txt' % D_s, 'w'))
    
    return demographic_events
