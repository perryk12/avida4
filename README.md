avida4
======

# How to run #
1. Figure out your executable name. For example, `logic9` or `mt`. 
2. `./mt` will run the executable
3. To set the config file, use:  `./exec -c config_file_name.cfg`. For example, `./mt -c major_transitions.cfg` will run the major transitions executable with the setting in the config file. 
4. To override one value on the command line, do something along the lines of the following; 
`./exec -c config_file_name.cfg --full.option.name=new_value`
For example, 
`./mt -c major_transitions.cfg --ea.gls.and_mutation_mult=6`
will set the value of `and_mutation_mult` to 6, overriding whatever was in the config file. 
5. To continue an existing run: 
`./exec -l <path_to_checkpoint_file>/checkpoint-1000000.xml.gz`
For example, `./mt -l /mnt/home/hjg/mt/033-gls-ramped/a_33/checkpoint-1000000.xml.gz`
6. To perform further analysis, load the check file and then run an analysis tool. 
`./mt_lr_gls -l /mnt/home/hjg/mt/033-gls-ramped/a_33/checkpoint-1000000.xml.gz --analyze lod_fitness --ea.analysis.input.filename /mnt/home/hjg/mt/033-gls-ramped/a_33/lod-1000000.xml.gz`
Note that this example, also loads a line of descent file which can also be analyzed. 

# Getting Started #
1. Try running an existing executable. 
2. Try changing a few values in the config file and rerunning. 
3. Figure out where the data you are looking at is being printed from. Hint: It's normally an event. In the case of the major transitions code, most of it is in the `mt_gls_propagule` which can be found in `mt.h`. 
