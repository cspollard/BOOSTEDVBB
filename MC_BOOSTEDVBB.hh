// -*- C++ -*-
#ifndef RIVET_MC_BOOSTEDVBB_HH
#define RIVET_MC_BOOSTEDVBB_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    typedef pair<string, double> JetCollection;

    class MC_BOOSTEDVBB : public Analysis {
        public:
            /// Constructor
            MC_BOOSTEDVBB()
                : Analysis("MC_BOOSTEDVBB") { }

            /// Book histograms and initialise projections before the run
            void init();

            /// Perform the per-event analysis
            void analyze(const Event& event);

            /// Normalise histograms etc., after the run
            void finalize();

        private:

            void bookChannel(const string& channel);

            void bookHistoAllChannels(const string& label,
                    const string& xlabel, int nxbins, double xmin, double xmax);

            void bookHistoAllChannels(const string& label,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel, int nybins, double ymin, double ymax);

            void bookProfileAllChannels(const string& label,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel);


            // map of channel, variable
            std::map<std::string, std::map<string, Histo1DPtr> >
                histo1DMap;

            std::map<std::string, std::map<string, Histo2DPtr> >
                histo2DMap;

            std::map<std::string, std::map<string, Profile1DPtr> >
                profile1DMap;

            std::vector<std::string> channels;
    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDVBB);
}

#endif
