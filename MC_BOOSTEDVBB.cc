// -*- C++ -*-
#include "MC_BOOSTEDVBB.hh"

#include <iostream>
#include <map>
#include <string>

#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/VariableR.hh"

using std::map;
using std::string;
using namespace Rivet::Cuts;

namespace Rivet {


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void MC_BOOSTEDVBB::init() {

        bookChannel("vvbb_inc");
        bookChannel("vvbb_0tag");
        bookChannel("vvbb_1tag");
        bookChannel("vvbb_2tag");

        bookChannel("lvbb_inc");
        bookChannel("lvbb_0tag");
        bookChannel("lvbb_1tag");
        bookChannel("lvbb_2tag");

        bookChannel("llbb_inc");
        bookChannel("llbb_0tag");
        bookChannel("llbb_1tag");
        bookChannel("llbb_2tag");

        bookHistoAllChannels("pTV", "V ${p_{\\rm T}}$ / GeV",
                50, 0, 1000*GeV);
        bookHistoAllChannels("pTJ", "jet ${p_{\\rm T}}$ / GeV",
                50, 0, 1000*GeV);
        bookHistoAllChannels("mJ", "jet mass / GeV",
                50, 0, 300*GeV);
        bookHistoAllChannels("dRbb", "$\\Delta R(b, b)$",
                50, 0, 4);


        FinalState allParticles(Cuts::abseta < 4.0);

        PromptFinalState promptLeptons = PromptFinalState(ChargedLeptons(allParticles));
        promptLeptons.acceptTauDecays(true);

        DressedLeptons leptonFinder(allParticles, promptLeptons,
                0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);

        addProjection(leptonFinder, "LeptonFinder");


        // calo jets constituents
        // INCLUDING NONPROMPT NEUTRINOS FOR NOW
        VetoedFinalState caloJetParts(allParticles);
        caloJetParts.addVetoOnThisFinalState(leptonFinder);
        FastJets fj10(caloJetParts, FastJets::ANTIKT, 1.0);
        fj10.useInvisibles(true);
        addProjection(fj10, "AKTCalo10");

        // track jets constituents
        ChargedFinalState trackParts(Cuts::abseta < 2.5 && Cuts::pT > 0.5*GeV);

        VetoedFinalState trackJetParts(trackParts);
        trackJetParts.addVetoOnThisFinalState(leptonFinder);
        addProjection(FastJets(trackJetParts, FastJets::ANTIKT, 0.2), "AKTTrack02");

        ZFinder zelelFinder(trackParts, Cuts::pT > 25*GeV,
                PID::ELECTRON, 71.2*GeV, 111.2*GeV);

        ZFinder zmumuFinder(trackParts, Cuts::pT > 25*GeV,
                PID::MUON, 71.2*GeV, 111.2*GeV);

        addProjection(zelelFinder, "ZelelFinder");
        addProjection(zmumuFinder, "ZmumuFinder");

        MissingMomentum missingMomentum(allParticles);
        addProjection(missingMomentum, "MissingMomentum");

        return;
    }


    /// Perform the per-event analysis
    void MC_BOOSTEDVBB::analyze(const Event& event) {

        const double w = event.weight();
        const vector<DressedLepton>& leptons =
            applyProjection<DressedLeptons>(event, "LeptonFinder").dressedLeptons();

        // TODO
        // does this have the correct sign?
        // update to Rivet2.4.0
        Vector3 met3V =
            applyProjection<MissingMomentum>(event, "MissingMomentum").vectorEt();
        met3V.setZ(0);

        FourMomentum met4V(met3V.mod2(), met3V.x(), met3V.y(), met3V.z());

        // find the appropriate nleptons channel.
        string channel;
        FourMomentum vboson;
        if (leptons.size() == 0) {
            channel = "vvbb";
            vboson = met4V;
        } else if (leptons.size() == 1) {
            channel = "lvbb";
            vboson = leptons.at(0).momentum() + met4V;
        } else if (leptons.size() == 2) {
            channel = "llbb";
            Particles zelels = applyProjection<ZFinder>(event, "ZelelFinder").bosons();
            Particles zmumus = applyProjection<ZFinder>(event, "ZmumuFinder").bosons();

            if (zelels.size() == 1 && zmumus.size() == 0)
                vboson = zelels.at(0);
            else if (zelels.size() == 0 && zmumus.size() == 1)
                vboson = zmumus.at(0);
            else
                vetoEvent;
        } else
            vetoEvent;

        const Jets& fatjets =
            applyProjection<FastJets>(event, "AKTCalo10").jetsByPt(Cuts::abseta < 2.0 && Cuts::pT > 250*GeV);

        // TODO
        // too harsh?
        // require exactly one high-pt large-R jet
        if (fatjets.size() != 1)
            vetoEvent;

        const Jet& fatjet = fatjets.at(0);

        const Jets& trackjets =
            applyProjection<FastJets>(event, "AKTTrack02").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 7.0*GeV);

        Jets matchedTrackJets;
        foreach (const Jet& trackjet, trackjets) {
            if (deltaR(fatjet, trackjet) < 1.0)
                matchedTrackJets.push_back(trackjet);
        }

        if (matchedTrackJets.size() < 2)
            vetoEvent;

        histo1DMap[channel + "_inc"]["pTV"]->fill(vboson.pT(), w);
        histo1DMap[channel + "_inc"]["pTJ"]->fill(fatjet.pT(), w);
        histo1DMap[channel + "_inc"]["mJ"]->fill(fatjet.mass(), w);
        histo1DMap[channel + "_inc"]["dRbb"]->fill(
                deltaR(matchedTrackJets.at(0), matchedTrackJets.at(1)), w);

        int nBtags = 0;
        foreach (const Jet& trackjet, matchedTrackJets) {
            // TODO
            // update for Rivet 2.4.0+
            // nBtags += (trackjet.bTagged(Cuts::pT > 5*GeV));
            nBtags += trackjet.bTagged();
        }

        if (nBtags == 0)
            channel += "_0tag";
        else if (nBtags == 1)
            channel += "_1tag";
        else
            channel += "_2tag";

        histo1DMap[channel]["pTV"]->fill(vboson.pT(), w);
        histo1DMap[channel]["pTJ"]->fill(fatjet.pT(), w);
        histo1DMap[channel]["mJ"]->fill(fatjet.mass(), w);
        histo1DMap[channel]["dRbb"]->fill(
                deltaR(matchedTrackJets.at(0), matchedTrackJets.at(1)), w);

        return;
    }


    /// Normalise histograms etc., after the run
    void MC_BOOSTEDVBB::finalize() {

        // normalize to 1/pb
        double norm = crossSection()/sumOfWeights();
        for (map< string, map<string, Histo1DPtr> >::iterator p = histo1DMap.begin(); p != histo1DMap.end(); ++p) {
            for (map<string, Histo1DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                q->second->scaleW(norm);
            }
        }

        for (map< string, map<string, Histo2DPtr> >::iterator p = histo2DMap.begin(); p != histo2DMap.end(); ++p) {
            for (map<string, Histo2DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                q->second->scaleW(norm);
            }
        }

        for (map< string, map<string, Profile1DPtr> >::iterator p = profile1DMap.begin(); p != profile1DMap.end(); ++p) {
            for (map<string, Profile1DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                q->second->scaleW(norm);
            }
        }

        return;
    }


    void MC_BOOSTEDVBB::bookChannel(const string& channel) {

        channels.push_back(channel);
        histo1DMap[channel] = map<string, Histo1DPtr>();
        histo2DMap[channel] = map<string, Histo2DPtr>();
        profile1DMap[channel] = map<string, Profile1DPtr>();

        return;
    }


    void MC_BOOSTEDVBB::bookHistoAllChannels(const string& label,
            const string& xlabel, int nxbins, double xmin, double xmax) {

        double xbinwidth = (xmax - xmin)/nxbins;

        char buff[100];
        sprintf(buff, "events / %.2f", xbinwidth);
        string ylabel = buff;

        foreach (const string& channel, channels) {
            histo1DMap[channel][label] =
                bookHisto1D(channel + "_" + label, nxbins, xmin, xmax, channel, xlabel, ylabel);
        }

        return;
    }


    void MC_BOOSTEDVBB::bookHistoAllChannels(const string& label,
            const string& xlabel, int nxbins, double xmin, double xmax,
            const string& ylabel, int nybins, double ymin, double ymax) {

        double xbinwidth = (xmax - xmin)/nxbins;
        double ybinwidth = (ymax - ymin)/nybins;

        char buff[100];
        sprintf(buff, "events / %.2f / %.2f", xbinwidth, ybinwidth);
        string zlabel = buff;

        foreach (const string& channel, channels) {
            histo2DMap[channel][label] =
                bookHisto2D(channel + "_" + label, nxbins, xmin, xmax,
                        nybins, ymin, ymax, channel,
                        xlabel, ylabel, zlabel);
        }

        return;
    }

    void MC_BOOSTEDVBB::bookProfileAllChannels(const string& label,
            const string& xlabel, int nxbins, double xmin, double xmax,
            const string& ylabel) {

        foreach (const string& channel, channels) {
            profile1DMap[channel][label] =
                bookProfile1D(channel + "_" + label, nxbins, xmin, xmax,
                        channel, xlabel, ylabel);
        }

        return;
    }

    //@}

} // Rivet
