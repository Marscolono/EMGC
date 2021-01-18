README

The following are details of the supplementary MATLAB code and data for the paper: Yuyanyuan Chen, Mingmin Xu, Xiaodan Fan, Cong Pian(2021),'Identifying functional modules using energy minimization with graph cuts'
  Contact: piancong@njau.edu.cn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMPILING THE SOURCE CODE

For energy minimisation we use the gco-v3.0 library available at http://vision.csd.uwo.ca/code/
To compile their source code use the following command in MATLAB: 
>> GCO_UnitTest
See their documentation for full details. 

Mex must be setup for both a C and C++ compiler. 
To check available compilers use the following command in MATLAB:
>> mex -setup
Note that for Mac OS X the only option for C and C++ is Xcode with Clang.  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USE THE CODE 

Genomic networks are represented as adjacency matrices. 

Whether or not the adjacency matrix is upper triangular is just a factor of 2 in beta. 
To find the value of beta^*, first overshoot and then trim the output of a second run as necessary. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

POSSIBLE ISSUES AND TROUBLESHOOTING 

Each of the .m files should return the same results as given in the paper. 
Make sure you can first reproduce the output before attempting to modify the code. 

The GCO code requires integers so scale the unary and pairwise potentials (DataCost and Neighbours arrays) as necessary. 

If the output is unexpected or cannot be produced, always just try something simpler first. 

Checking which label is dominant can be done using the following code: 

%---------------------------------------------------
num_hits = NaN*ones(nLabels,size(min_energy_labels,2));

for i = 1:size(min_energy_labels,2)
    for l = 1:nLabels
        num_hits(l,i) = sum(min_energy_labels(:,i) == l);
    end
end

plot(num_hits')
%---------------------------------------------------

Note that isolated nodes will never have the dominant label unless they 
started with it and so take that into account when determining dominance. 