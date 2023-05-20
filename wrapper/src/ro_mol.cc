#include "rust/cxx.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

namespace RDKit {
    using ExplicitBitVect = ::ExplicitBitVect;

    std::shared_ptr<ROMol> copy_mol(std::shared_ptr<ROMol> mol) {
        return std::shared_ptr<ROMol>(new ROMol(*mol));
    }

    std::shared_ptr<ROMol> smiles_to_mol(const std::string &smiles) {
        ROMol *mol = SmilesToMol(smiles);

        return std::shared_ptr<ROMol>(mol);
    }

    std::shared_ptr<ROMol> smiles_to_mol_with_params(const std::string &smiles, std::shared_ptr<SmilesParserParams> params) {
        ROMol *mol = SmilesToMol(smiles, *params);

        return std::shared_ptr<ROMol>(mol);
    }
    std::shared_ptr<SmilesParserParams> new_smiles_parser_params() {
        return std::shared_ptr<SmilesParserParams>(new SmilesParserParams());
    }
    void smiles_parser_params_set_sanitize(std::shared_ptr<SmilesParserParams> params, bool sanitize) {
        params->sanitize = sanitize;
    }

    rust::String mol_to_smiles(std::shared_ptr<ROMol> mol) {
        return MolToSmiles(*mol);
    }

    std::unique_ptr<std::vector<std::string>> detect_chemistry_problems(std::shared_ptr<ROMol> mol) {
        std::vector<std::unique_ptr<RDKit::MolSanitizeException>> exceptions = RDKit::MolOps::detectChemistryProblems(*mol);
        std::vector<std::string> *get_types = new std::vector<std::string>;
        get_types->reserve(exceptions.size());

        for (auto &t: exceptions) {
            get_types->push_back(t->getType());
        }

        return std::unique_ptr<std::vector<std::string>>(get_types);
    }

    rust::String draw_mol(std::shared_ptr<ROMol> mol) {
        RDKit::MolDraw2DSVG drawer(500, 400, -1, -1, true);
        drawer.drawOptions().additionalAtomLabelPadding = 0.05;
        drawer.drawOptions().explicitMethyl = true;
        drawer.drawOptions().padding = 0.0;
        MolDraw2DUtils::prepareAndDrawMolecule(drawer, *mol);
        drawer.finishDrawing();
        rust::String data = drawer.getDrawingText();
        return data;
    }
}