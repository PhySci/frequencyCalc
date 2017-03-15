The MatLab class "freqCalc" provides functionality for calculation of local FMR frequencies, using a semi-analitical approach described by C. S. Davies, V. D. Poimanov, and V. V. Kruglyak (publication in progress).

## How to use it?

1. Clone the repository (or download and unpack) in a folder.
2. Add the folder to MatLab search path.
3. Calculate a ground state for the system under study, save it as a "static.stc"
4. Run MatLab, go to a folder with the ground state file "static.stc".
5. Create instance of freqCalc class (here and below I'm using name "FC" for the instance)


    FC = freqCalc;
	
6. Calculate deflected states of static magnetisation


    FC.calcDeflection();
 	
7. Now you have to calculate magnetic fields for the given states using a micromagnetic package. If you use muMax3, then open a file with micromagnetic model (a description of the system), remove relaxation and any excitations (CW, pulse etc) and add next lines:
    
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
	

First four lines load the ground state "static.stc" and save demagnetization, exchange and effective fields. Next four lines load the deflectated state "outplane.stc" and save fields. Last four lines load the deflectated state "inplance.stc" and save fields.

8. Run the simulation. It takes a few minutes. As a results you will have a plenty of files of magnetic fields.

9. Go to MatLab. Set correct value of bias magnetic field (for instance, H = 500 Oe).
      FC.H = 500;

10. Run calculatation of demagnetization tensor and spatial map of FMR frequency:

     FC.calcDemagTensor('field','demag');

There are a few options which fields will be taken into account:
    
     FC.calcDemagTensor('field','demag') - calculate demagnetization tensor using demagnetization field only (default);
    
     FC.calcDemagTensor('field','demag+exchange') - calculate demagnetization tensor using demagnetization and exchnage fields only (default); 

Output of the script is a few images (diagonal componens of the demagnetization tensor, map of frequency, local fields). But you can read all these values as properies of the class instance (for example, FC.Nxx is three dimensional array of Nxx demagnetization component, FC.freq is three-dimensional array of FMR frequency and so on). 
	
