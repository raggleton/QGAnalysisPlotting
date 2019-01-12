import math

for boson,name in [("w","enj"),("z","eej")]:
  f=open(name+".dat")
  pt_bins=[]
  values={}
  entry=None
  for l in f.readlines():
    split=l.strip().split(" ")
    if name+"_pTV" in l:
       entry=split[-1]
       values[entry]=[]
    if "#" in l or entry==None or len(split)<5: continue
    if entry==name+"_pTV_LO":
       pt_bins+=[float(split[0])]
    values[entry]+=[float(split[2])]

  print(pt_bins)
  kNNLO=[values[name+"_pTV_NNLO"][i]/values[name+"_pTV_LO"][i] for i in range(len(pt_bins))]
  print(kNNLO)
  d1kNLO=[values[name+"_pTV_d1K_NLO"][i] for i in range(len(pt_bins))]
  print(d1kNLO)
  d2kNLO=[values[name+"_pTV_d2K_NLO"][i] for i in range(len(pt_bins))]
  print(d2kNLO)
  d1kappaEW=[values[name+"_pTV_d1kappa_EW"][i] for i in range(len(pt_bins))]
  print(d1kappaEW)
  d2kappaEW=[values[name+"_pTV_d2kappa_EW"][i] for i in range(len(pt_bins))]
  print(d2kappaEW)
  dPDF=[]
  for i in range(len(pt_bins)):
     mean=0
     count=0
     for pdfname in values.keys():
       if "PDF" in pdfname:
          mean+=values[pdfname][i]
          count+=1
     mean/=count
     rms=0
     for pdfname in values.keys():
       if "PDF" in pdfname:
          rms+=pow(values[pdfname][i]-mean,2)
     rms=math.sqrt(rms/count)
     dPDF+=[rms]
  print(dPDF)

  import ROOT
  import array

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",kNNLO)
  gkNNLO=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gkNNLO.SetName("kNNLO")

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",d1kNLO)
  gd1kNLO=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gd1kNLO.SetName("d1kNLO")

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",d2kNLO)
  gd2kNLO=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gd2kNLO.SetName("d2kNLO")

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",d1kappaEW)
  gd1kappaEW=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gd1kappaEW.SetName("d1kappaEW")

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",d2kappaEW)
  gd2kappaEW=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gd2kappaEW.SetName("d2kappaEW")

  pointsX=array.array("f",pt_bins)
  pointsY=array.array("f",dPDF)
  gdPDF=ROOT.TGraph(len(pointsX),pointsX,pointsY)
  gdPDF.SetName("dPDF")

  out=ROOT.TFile(boson+"_pt_spectrum.root","RECREATE")
  gkNNLO.Write()
  gd1kNLO.Write()
  gd2kNLO.Write()
  gd1kappaEW.Write()
  gd2kappaEW.Write()
  gdPDF.Write()
  out.Close()
  