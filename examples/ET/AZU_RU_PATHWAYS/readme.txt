1. Start HARLEM and Load PDB file of the ET complex 
( Azurin Ruthenated at HIS 83 ): 

   a) Main Menu: File->Open Molecule
   b) Choose Directory: c:\harlem\examples\ET\AZU_RU_PATHWAYS 
   c) File Type:  Prtein Data Bank
   d) Choose and Load File azu_83_ru.pdb

You can also drop the icon of PDB file azu_83_ru.pdb on harlem.exe icon. 

2. To see metal centers easier restrict view to atoms within 4.0 And from 
   residues CU ( copper center) and RBP (ruthenium bipyridine complex) 
   
   a) In command prompt type:  restrict within(4.0,CU) or within(4.0,RBP) 
   b) press "Execute COmmand"
   c) Display->Center View Sel Atoms   - to expand view to selected atoms

3. Specify donor and acceptor in Edit Redox Centers Dialog:
 
a) Open dialog by Application->ET->'Edit Redox Centers' 
b) Choose DONOR toggle button
c) Click on Copper atom 
    In Expr Window appear:  $AZU_83_RU$CU129:A.CU
d) Add Copper atom to the DONOR group  pressing "Add From Expr" 
e) Choose ACCEPTOR toggle button
f) Click on Copper atom 
    In Expr Window appear:  $AZU_83_RU$RBP130:A.RU
g) Add Ruthenium atom to the ACCEPTOR group  pressing "Add From Expr"
h) Toggle between DONOR and ACCEPTOR groups:
    Atoms of the groups appear as  small spheres. 
    All atoms of the donor(acceptor) assumed coupled with the value 1.0 (short-circuited).

4. Turn back to All Atom display: 
    a) Main Menu:  Edit->Select All 
    b) Display->Display Mode->WireFrame 
    
5. Run Pathways:

    a) Applications->ET->"Run Pathways" to open PATHWAYS Dialog.
    b)  One can do three types of PATHWAYS calculations:
    i) Best Path calculations

       Press "Compute Best Path" 
       Atom Composition of the BEST PATH between specified DONOR and ACCEPTOR groups 
       will appear in HARLEM CONSOLE Window. Best Path will be seen as yellow sticks on the molecular view.

    ii) Coupling Map: Coupling of all atoms to the donor

       "Compute ET Coupling Map"  Coupling map will be printed to HARLEM Console Window. and 
	Protein will be coloured by coupling values (red- max, blue- min)

        One can save the molecule in pdb format: File->Save Molecule 
								 Choose File Type : Protein Data Bank
								Choose File Name and press "Save File"

	The computed log of coupling values will be in the temperature factors of the atoms in the pdb file.

   iii) 3-rd type of PATHWAYS calculations is the selection of all atoms which belong 
	to the paths which have coupling values within certain threshold value from the BEST PATH 
	say within a factor of 0.5 or 0.1. These calculations give you an estimate which residues 
	are important for mediating ET coupling between Donor and Acceptor.

        1) Choose: Threshol coupling relative to the best path
	2) press "Select ET Coupled Atoms"
           Atoms on the important paths and donor and acceptor will be plotted as tubes. 
	   One can get idea what atoms are important by clicking on selected atoms - their identifiers 
	   will appear in HARLEM CONSOLE window.
	   
	   