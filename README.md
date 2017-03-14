The MatLab class "freqCalc" provides functionality for calculation of local FMR frequencies, using a semi-analitical approach described by C. S. Davies, V. D. Poimanov, and V. V. Kruglyak (publication in progress).

How to use it?

1) Clone the repository (or download and unpack) in a folder.
2) Add the folder to MatLab search path.
3) Calculate a ground state for the system under study, save it as a "static.stc"
4) Run MatLab, go to a folder with the ground state file "static.stc".
5) Create instance of freqCalc class (here and below I'm using name "FC" for the instance) 
 

  

// static
m.loadFile("static.stc");  
saveAs(B_demag,"B_demag_static");
saveAs(B_exch,"B_exch_static"); 
saveAs(B_eff,"B_eff_static"); 
 
// out-of-plane deflection
m.loadFile("outplane.stc");
saveAs(B_demag,"B_demag_outplane");
saveAs(B_exch,"B_exch_outplane");
saveAs(B_eff,"B_eff_outplane");

// in-plane deflection 
m.loadFile("inplane.stc");  
saveAs(B_demag,"B_demag_inplane");
saveAs(B_exch,"B_exch_inplane");
saveAs(B_eff,"B_eff_inplane");