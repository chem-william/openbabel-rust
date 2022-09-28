#pragma once
#include <memory>
#include "rust/cxx.h"
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>

namespace OpenBabel {
    class OBMol;
    class OBSmartsPattern;
    class OBConversion;
    class OBForceField;

    // Debug
    void print_global_instances();

    // OBConversion
    // std::unique_ptr<OBMol> OBConversion_smi_to_mol(const std::string &smiles);

    // OBForceField
    OBForceField* OBForceField_find_forcefield(const std::string &ff_name);
    unsigned int OBForceField_setup(const std::unique_ptr<OBMol> & pMol, OBForceField* pFF);
    void OBForceField_conjugate_gradients(OBForceField* pFF, u_int32_t steps, double econv);
    void OBForceField_conjugate_gradients_initialize(OBForceField* pFF, u_int32_t steps, double econv);
    bool OBForceField_conjugate_gradients_take_n_steps(OBForceField* pFF, u_int32_t n);

    void OBForceField_steepest_descent(OBForceField* pFF, u_int32_t steps, double econv);
    void OBForceField_steepest_descent_initialize(OBForceField* pFF, u_int32_t steps, double econv);
    bool OBForceField_steepest_descent_take_n_steps(OBForceField* pFF, u_int32_t n);
    double OBForceField_energy(OBForceField* pFF);


    // OBMol
    std::unique_ptr<OBMol> OBMol_from_smiles(const std::string &smiles);
    unsigned int OBMol_num_atoms(const std::unique_ptr<OBMol> & pMol);
    unsigned int OBMol_num_bonds(const std::unique_ptr<OBMol> & pMol);
    unsigned int OBMol_num_hvy_atoms(const std::unique_ptr<OBMol> & pMol);
    double OBMol_get_mol_wt(const std::unique_ptr<OBMol> & pMol);

    // OBFingerprint
    typedef std::vector<unsigned int> FPData;
    std::unique_ptr<FPData> OBFingerprint_get_fingerprint(const std::string &fp_name, const std::unique_ptr<OBMol> & pMol, u_int32_t nbits);
    // std::unique_ptr<FPData> OBFingerprint_get_fingerprint_in_batch(const std::string &fp_thread_name, const rust::Vec<rust::String> & smiles_vec, u_int32_t nbits);
    // deprecated: slow performance, root cause to be identified

    // OBSmartsPattern
    std::unique_ptr<OBSmartsPattern> OBSmartsPattern_from_smarts(const std::string &smarts);
    unsigned int OBSmartsPattern_num_atoms(const std::unique_ptr<OBSmartsPattern> & pSP);
    unsigned int OBSmartsPattern_num_bonds(const std::unique_ptr<OBSmartsPattern> & pSP);
    unsigned int OBSmartsPattern_num_matches(const std::unique_ptr<OBSmartsPattern> & pSP);
    std::unique_ptr<std::vector<int>> OBSmartsPattern_match(const std::unique_ptr<OBSmartsPattern> & pSP, const std::unique_ptr<OBMol> &pMol);
}