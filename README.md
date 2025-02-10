# OrderedWSets
The source code written in C++ for Ordered WSets.

The source code written in R for processing datasets.

FSMs are given in the following format:

FSM_ID #_States #_Transitions #_Inputs #_Outputs Void

initial_state destination_state input output

...........................................

FSMsStat and benchFSMsStat files have the following columns:

FSM_ID: FSM identifier.

FSMStates: Number of states.

FSMInputs: Number of inputs.

FSMOutputs: Number of outputs.

FSMSize: Number of transitions.

W_Time: Time required to construct W set.

W_Memory: Memory required to construct W set.

SIS_Tr: Total number of inputs used for transfer in the state identification.

SIS_Size: Total number of sequences in the state identifier.

SIS_Cost: Total number of inputs of the state identifier.

TS_Size: Total number of sequences of the test suite.

TS_Cost: The total number of inputs of the test suite. 

Algorithm: Algorithm name.


FileBindings file maps the names of the benchmark FSMs given in fileBindings.txt with the FSM id's given in benchFSMs.txt. Please note in BenchFSMs.txt every FSM has an FSM_ID. So FSM_ID = 0 is for bbara_minimized.dot.


To compile the Source Code, create a project under Visual Studio 19 and store the source and txt files. If you fail to compile the code, please contact Uraz via u.turker@lancaster.ac.uk.
