var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Epidemics","category":"page"},{"location":"#Epidemics","page":"Home","title":"Epidemics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Epidemics.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Epidemics]","category":"page"},{"location":"#Epidemics.Intervention","page":"Home","title":"Epidemics.Intervention","text":"Npi(timebegin, timeend, contact_reduction)\n\nA structure to hold the end points and strength of a non-pharmaceutical intervention.\n\n\n\n\n\n","category":"type"},{"location":"#Epidemics.Population","page":"Home","title":"Epidemics.Population","text":"Population(name, demographyvector, initialconditions, contact_matrix)\n\nA structure to hold population characteristics, including:\n\n'name': A name for the population.\n'demography_vector': A numeric vector of the number of individuals in each\n\nage or demographic group of the population.\n\n'initial_conditions': A numeric matrix representing the proportions of each\n\nage or demographic group that are in one of the epidemiological compartments.\n\n'contact_matrix': A matrix giving the contacts between the demographic groups\n\nin the population. Must be a square matrix.\n\n\n\n\n\n","category":"type"},{"location":"#Epidemics.epidemic-Tuple{}","page":"Home","title":"Epidemics.epidemic","text":"epidemic(model_name, init, contact_matrix, demography_vector, r0,\n    preinfectious_period, infectious_period, interv, time_end, increment\n)\n\nModel the progression of an epidemic, with age- or demographic-group specific contact patterns and proportions, epidemiological parameters, and interventions.\n\n\n\n\n\n","category":"method"},{"location":"#Epidemics.seir!-NTuple{4, Any}","page":"Home","title":"Epidemics.seir!","text":"seir!(du, u, parameters, t)\n\nA simple SEIR epidemic model function that allows for multiple demographic     groups. This function is intended to be called internally from     epidemic.     The function expects the parameters argument to be a four element vector     with the following elements:     - a vector of beta, the transmission rate, where each element of the     vector represents the transmission rate of the pathogen within a specific     demographic group;     - a vector of alpha, the group-specific rate of conversion from exposed      to infectious;     - a vector of gamma, the group-specific rate of recovery;     - a matrix specifying the contacts between demographic groups;     - a matrix of the interventio applied to each age group, see Npi\n\n\n\n\n\n","category":"method"}]
}
