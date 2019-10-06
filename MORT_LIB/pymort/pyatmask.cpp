    atomvec_t (*ambmsk_mask_atom_fn)(const molecule_t&, const string& ) = &mask_atom;
    atomvec_t (*smarts_mask_atom_fn)(const molecule_t&, const string& ) = &smarts_mask_atom;
    def("_ambmsk_mask_atom", ambmsk_mask_atom_fn);
    def("_smarts_mask_atom", smarts_mask_atom_fn);
 
