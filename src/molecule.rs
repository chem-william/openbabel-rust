//! Openbabel Molecule
//! 
//! 
//! # Create a Molecule from a SMILES string 
//! ```
//! use openbabel::molecule;
//! 
//! let mol = molecule::Molecule::new_from_smiles("c1ccccc1");
//! assert_eq!(mol.num_atoms(), 6);
//! ```


use ob_rs::ob;

pub struct Molecule {
    pub ob_mol: cxx::UniquePtr<ob::OBMol>
}

impl Molecule {
    pub fn new_from_smiles(smiles: &str) -> Self {
        cxx::let_cxx_string!(smiles_cxx = smiles);
        Self { ob_mol: ob::OBMol_from_smiles(&smiles_cxx) }
    }
    
    pub fn is_valid(&self) -> bool { !self.ob_mol.is_null() }
    pub fn num_atoms(&self) -> u32 { ob::OBMol_num_atoms(&self.ob_mol) }
    pub fn num_bonds(&self) -> u32 { ob::OBMol_num_bonds(&self.ob_mol) }
    pub fn num_hvy_atoms(&self) -> u32 { ob::OBMol_num_hvy_atoms(&self.ob_mol) }
    pub fn get_mol_wt(&self) -> f64 { ob::OBMol_get_mol_wt(&self.ob_mol) }
    pub fn num_rotors(&self) -> u32 { ob::OBMol_num_rotors(&self.ob_mol) }
    pub fn get_formula(&self) -> String { ob::OBMol_get_formula(&self.ob_mol) }
    pub fn get_coordinates(&self) -> Vec<[f64; 3]> {
        ob::OBMol_get_coordinates(&self.ob_mol).chunks_exact(3).map(|v| [v[0], v[1], v[2]]).collect()
    }
}

#[cfg(test)]
mod test_mod_molecule {
    use super::*;

    #[test]
    fn test_molecule() {
        let mol = Molecule::new_from_smiles("CCNCC");
        assert!(mol.is_valid());
        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 4);
        assert_eq!(mol.num_hvy_atoms(), 5);
        assert_eq!(mol.num_rotors(), 2);
        assert_eq!(mol.get_formula(), "C4H11N");

        assert!((mol.get_mol_wt() - 73.137).abs() < 1e-3, "mol wt is {}", mol.get_mol_wt());
    }
}
