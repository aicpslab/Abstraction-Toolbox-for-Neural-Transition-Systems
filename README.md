# Abstraction-Toolbox-for-Neural-Transition-Systems
This is the code for the abstraction of neural transition systems for the learning model of the maglev model. The description for the maglev model can be found at https://ww2.mathworks.cn/help/deeplearning/ug/maglev-modeling.html, while the neural network approximation of its dynamics is at https://github.com/xiangweiming/ignnv.

# Related tools and software
This toolbox makes use of IGNNV (at:https://github.com/xiangweiming/ignnv) for reachability analysis; UPPAAL 5.0 (at: https://uppaal.org/)) for LTL specification verification.

# Setup
1. MATLAB pre-installed is required.
2. UPPAAL 5.0 pre-install is required.
3. Add the folder "Abstraction-Toolbox-for-Neural-Transition-Systems" to the MATLAB workpath.

# Running tests and examples

1. Open .../Abstraction-Toolbox-for-Neural-Transition-Systems/Maglev_example.mlx with MATLAB to execute the scripts for testing/analyzing examples and generate the model graph named "MaglevModelGraph1.xml".
2. Open .../MaglevModelGraph1.xml with UPPAAL to verify if the neural transition system satisfy certain specifications, and return counterexample if not.
