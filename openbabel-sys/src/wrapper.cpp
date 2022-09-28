#include <sstream>
#include <vector>
#include <openbabel/fingerprint.h>
#include <openbabel/oberror.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include "wrapper.h"

namespace OpenBabel {

// OBConversion 

// unique_ptr<OBMol> OBConversion_smi_to_mol(const string &smiles) {
//     OBSmilesParser ob_sp = OBSmilesParser();
//     OBMol mol = OBMol();
//     if (ob_sp.SmiToMol(mol, smiles)) {
//         return make_unique<OBMol>(move(mol));
//     } else {
//         return unique_ptr<OBMol>(nullptr);
//     }
// }

unsigned int OBConversion_set_in_format(const unique_ptr<OBConversion> & pConv, const string &input_format) {
    OBFormat* pFormat = OBConversion::FindFormat(input_format);

    if (!pConv.get()->SetInFormat(pFormat)) {
	stringstream errorMsg;
        errorMsg << "OBConversion::SetInFormat(" << input_format << ")" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return 1;
    }

    return 0;
}

unsigned int OBConversion_set_out_format(const unique_ptr<OBConversion> & pConv, const string &output_format) {
    OBFormat* pFormat = OBConversion::FindFormat(output_format);

    if (!pConv.get()->SetOutFormat(pFormat)) {
	stringstream errorMsg;
        errorMsg << "OBConversion::SetOutFormat(" << output_format << ")" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return 1;
    }

    return 0;
}

unsigned int OBConversion_set_in_and_out_formats(const unique_ptr<OBConversion> & pConv, const string &input_format, const string &output_format) {
    OBFormat* pInFormat = OBConversion::FindFormat(input_format);
    OBFormat* pOutFormat = OBConversion::FindFormat(output_format);

    if (!pConv.get()->SetInAndOutFormats(pInFormat, pOutFormat)) {
	stringstream errorMsg;
        errorMsg << "OBConversion::SetInAndOutFormats(" << output_format << ")" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return 1;
    }

    return 0;
}

unsigned int OBConversion_read_string(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol, const string &input) {
    if (!pConv.get()->ReadString(pMol.get(), input.c_str())) {
        stringstream errorMsg;
        errorMsg << "OBConversion::ReadString error" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
	return 1;
    }
    return 0;
}

rust::String OBConversion_write_string(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol) {
    return pConv.get()->WriteString(pMol.get());
}

unsigned int OBConversion_read_file(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol, const string &input_path) {
    if (!pConv.get()->ReadFile(pMol.get(), input_path)) {
        stringstream errorMsg;
        errorMsg << "OBConversion::ReadFile error" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
	return 1;
    }
    return 0;
}

rust::Vec<rust::String> OBConversion_get_supported_input_format() {
    OBConversion obconv = OBConversion();
    rust::Vec<rust::String> result {};
    for (auto element : obconv.GetSupportedInputFormat()) {
	result.push_back(element);
    }
    return result;
}

rust::Vec<rust::String> OBConversion_get_supported_output_format() {
    OBConversion obconv = OBConversion();
    rust::Vec<rust::String> result {};
    for (auto element : obconv.GetSupportedOutputFormat()) {
	result.push_back(element);
    }
    return result;
}

// OBConversion - End

// OBForceField
OBForceField* OBForceField_find_forcefield(const string &ff_name) {
    OBForceField* pFF = OBForceField::FindForceField(ff_name.c_str());

    if (!pFF) {
        stringstream errorMsg;
	errorMsg << "OBForceField::FindForceField error" << endl;
	obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }

    return pFF;
}

unsigned int OBForceField_setup(const unique_ptr<OBMol> & pMol, OBForceField* pFF) {
    if (!pFF->Setup(*pMol)) {
        stringstream errorMsg;
	errorMsg << "OBForceField->Setup() error" << endl;
	obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
	return 1;
    }
    return 0;
}

void OBForceField_conjugate_gradients(OBForceField* pFF, u_int32_t steps, double econv) {
    pFF->ConjugateGradients(steps, econv);
}

void OBForceField_conjugate_gradients_initialize(OBForceField* pFF, u_int32_t steps, double econv) {
    pFF->ConjugateGradientsInitialize(steps, econv);
}

bool OBForceField_conjugate_gradients_take_n_steps(OBForceField* pFF, u_int32_t n) {
    return pFF->ConjugateGradientsTakeNSteps(n);
}

void OBForceField_steepest_descent(OBForceField* pFF, u_int32_t steps, double econv) {
    pFF->SteepestDescent(steps, econv);
}

void OBForceField_steepest_descent_initialize(OBForceField* pFF, u_int32_t steps, double econv) {
    pFF->SteepestDescentInitialize(steps, econv);
}

bool OBForceField_steepest_descent_take_n_steps(OBForceField* pFF, u_int32_t n) {
    return pFF->SteepestDescentTakeNSteps(n);
}

double OBForceField_energy(OBForceField* pFF) { return pFF->Energy(); }
const string OBForceField_get_unit(OBForceField* pFF) { return pFF->GetUnit(); }

// OBForceField End


// OBMol
unique_ptr<OBMol> OBMol_new() { return unique_ptr<OBMol>(new OBMol()); }
unique_ptr<OBConversion> OBConversion_new() { return unique_ptr<OBConversion>(new OBConversion()); }
unique_ptr<OBMol> OBMol_from_smiles(const string &smiles) {
    unique_ptr<OBMol> pMol(new OBMol());
    stringstream ss(smiles);
    OBConversion conv(&ss);

    if (!conv.SetInFormat("smi")) {
        stringstream errorMsg;
        errorMsg << "OBConversion::SetInFormat (\"smi\")  error" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    }

    if(!conv.Read(pMol.get())) {
        stringstream errorMsg;
        errorMsg << "OBConversion::Read error" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return unique_ptr<OBMol>(nullptr);
    } 

    return pMol;
}

unsigned int OBMol_num_atoms(const unique_ptr<OBMol> & pMol) { return pMol->NumAtoms(); }
unsigned int OBMol_num_bonds(const unique_ptr<OBMol> & pMol) { return pMol->NumBonds(); }
unsigned int OBMol_num_hvy_atoms(const unique_ptr<OBMol> & pMol) { return pMol->NumHvyAtoms(); }
double OBMol_get_mol_wt(const unique_ptr<OBMol> & pMol) { return pMol->GetMolWt(); }

// OBMol End

// OBFingerprint

unique_ptr<FPData> OBFingerprint_get_fingerprint(const string &fp_name, const unique_ptr<OBMol> & pMol, u_int32_t nbits) {
    FPData fps;
    fps.resize(nbits / 32);
    OBFingerprint* pFP = OBFingerprint::FindFingerprint(fp_name.c_str());

    if (!pFP) {
        stringstream errorMsg;
        errorMsg << "Cannot find fingerprint " << fp_name << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        fill(fps.begin(), fps.end(), 0);
    } else {
        if(!pFP->GetFingerprint(pMol.get(), fps, nbits)) {
            stringstream errorMsg;
            errorMsg << "Error on generating fingerprint " << fp_name << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
            fill(fps.begin(), fps.end(), 0);
        }
    }

    return make_unique<FPData>(move(fps));
}

// OBFingerprint - End


// OBSmartsPattern

unique_ptr<OBSmartsPattern> OBSmartsPattern_from_smarts(const string &smarts) {
    unique_ptr<OBSmartsPattern> pSP(new OBSmartsPattern());
    pSP->Init(smarts);
    return pSP;
}

unsigned int OBSmartsPattern_num_atoms(const unique_ptr<OBSmartsPattern> & pSP) { return pSP->NumAtoms(); }
unsigned int OBSmartsPattern_num_bonds(const unique_ptr<OBSmartsPattern> & pSP) { return pSP->NumBonds(); }
unsigned int OBSmartsPattern_num_matches(const unique_ptr<OBSmartsPattern> & pSP) { return pSP->NumMatches(); }

unique_ptr<vector<int>> OBSmartsPattern_match(const unique_ptr<OBSmartsPattern> & pSP, const unique_ptr<OBMol> & pMol) {
    pSP->Match(*pMol);
    // CxxVector does not support nested C++ vector (vector<vector>)
    vector<int> result {};
    for (vector<vector<int>>::iterator i = pSP->GetMapList().begin(); i != pSP->GetMapList().end(); ++i) {
        result.insert(result.end(), i->begin(), i->end());
    }

    return make_unique<vector<int>>(move(result));
}

// OBSmartsPattern - End


} // namespace OpenBabel
