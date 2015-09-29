from AthenaCommon.AppMgr import ServiceMgr as svcMgr
theApp.EvtMax = -1

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = ["EVNT.05915947._001376.pool.root.1"]

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
rivet.Analyses += [ 'MC_BOOSTEDVBB' ]
rivet.RunName = ""
rivet.HistoFile = "boostedvbb.yoda"
import os
rivet.AnalysisPath = os.environ['PWD']
rivet.CrossSection = 1 
job += rivet 

from GaudiSvc.GaudiSvcConf import THistSvc

svcMgr += THistSvc() 
svcMgr.THistSvc.Output = ["Rivet DATAFILE='boostedvbb.root' OPT='RECREATE'"] 

job += rivet
