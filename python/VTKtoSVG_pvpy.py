import sys
#sys.path.append("C:\Program Files\ParaView 5.10.1-Windows-Python3.9-msvc2017-AMD64\\bin\Lib\site-packages")
#sys.path.append("C:\Program Files\ParaView 5.10.1-Windows-Python3.9-msvc2017-AMD64\\bin")

from paraview.simple import *

import time
import io
import os
import glob

files=glob.glob("../output/SVG files/Microenvironment/*.svg")
for f in files:
    os.remove(f)

Disconnect()
Connect()

files_vtr=glob.glob('../output/VTK files/microenvironment_*.vtr')
files_vtr.sort()

#print(files_vtr)

#RenderView1=GetActiveView()
#Delete(renderView1)

for x in range(0,len(files_vtr)):
    sx=str(x)
    if x<100:  
        sx="0"+sx
    if x<10:
        sx="0"+sx

    output_debris_path = '../output/SVG files/Microenvironment/Debris_' + sx + '.svg'
    output_TNF_path = '../output/SVG files/Microenvironment/TNFa_' + sx + '.svg'
    output_TGF_path = '../output/SVG files/Microenvironment/TGFb_' + sx + '.svg'
    output_IL10_path = '../output/SVG files/Microenvironment/IL10_' + sx + '.svg'
    output_IFNg_path = '../output/SVG files/Microenvironment/IFNg_' + sx + '.svg'
    input_path = '../output/VTK files/microenvironment_' + sx + '.vtr'
    directory = os.path.dirname(output_debris_path)

    # try:
        # os.makedirs(directory)
    # except OSError:
        # if not os.path.isdir(directory):
            # raise 

    # Open the corresponding VTK

    reader = OpenDataFile(input_path)

    RenderView1=CreateRenderView()
    RenderView1.UseColorPaletteForBackground=0
    SetViewProperties(Background = [1, 1, 1])
    RenderView1.OrientationAxesVisibility = 0

    Display1=Show()
    camera = GetActiveCamera()

    ##DEBRIS

    ColorBy(Display1,'Debris')
    Display1.RescaleTransferFunctionToDataRange(True)

    dataDEBColorMap = GetColorTransferFunction('Debris')
    dataDEBOpacityMap = GetOpacityTransferFunction('Debris')

    dataDEBColorMap.EnableOpacityMapping = True

    dataDEBColorMap.ApplyPreset('bone_Matlab', True)
    dataDEBColorMap.InvertTransferFunction()
    
    dataDEBColorMap.RescaleTransferFunction(0.0, 1e-4)  
    dataDEBOpacityMap.RescaleTransferFunction(0.0, 1e-4)

    Display1.SetScalarBarVisibility(RenderView1, True)
    DEBscalarBar=GetScalarBar(dataDEBColorMap,RenderView1)
    DEBscalarBar.TitleFontSize = 15
    DEBscalarBar.TitleColor=[0,0,0]
    DEBscalarBar.LabelFontSize=8 
    DEBscalarBar.LabelColor=[0,0,0]

    #RenderView1.UseColorPaletteForBackground="WhiteBackground"

    #DEBscalarBar.DrawSubTickMarks=False
    # TNFscalarBar.MaximumNumberofLabels = 3

    RenderView1.ResetCamera()

    # RenderView1 = GetRenderView()
    # RenderView1.CameraViewUp = [0, 0, 0]
    # RenderView1.CameraPosition = [0, 0, 640]
    # RenderView1.CameraClippingRange = [0, 0]

    camera = GetActiveCamera()
    #camera.Dolly(1.4)

    # camera.SetFocalPoint(0, 0, 0)
    # camera.SetPosition(0, 0, -10)
    # camera.SetViewUp(0, 1, 0)
    # camera.SetViewAngle(30)

    # renderView1 = GetActiveViewOrCreate('RenderView')
    # renderView1.ViewSize = [1670, 1091]

    # outDisplay = GetDisplayProperties(reader, view=renderView1)

    Render()
    
    ExportView(output_debris_path, view=RenderView1, Plottitle='ParaView GL2PS Export')

    Delete(DEBscalarBar)
    Delete(RenderView1)

    RenderView1 = CreateRenderView()
    RenderView1.UseColorPaletteForBackground=0
    SetViewProperties(Background = [1, 1, 1])
    RenderView1.OrientationAxesVisibility = 0

    Display1=Show()
    camera = GetActiveCamera()

    ##TNFa

    ColorBy(Display1,'TNFa')
    Display1.RescaleTransferFunctionToDataRange(True)

    # TNFdata = GetColorTransferFunction('debris')
    # TNFdata.ApplyPreset('Cold and Hot', True)

    dataTNFColorMap = GetColorTransferFunction('TNFa')
    dataTNFOpacityMap = GetOpacityTransferFunction('TNFa')

    dataTNFColorMap.EnableOpacityMapping = True

    dataTNFColorMap.ApplyPreset('Reds', True)
    dataTNFColorMap.InvertTransferFunction()
    
    dataTNFColorMap.RescaleTransferFunction(0.0, 2e-11)
    dataTNFOpacityMap.RescaleTransferFunction(0.0, 2e-11)

    Display1.SetScalarBarVisibility(RenderView1, True)
    TNFscalarBar=GetScalarBar(dataTNFColorMap,RenderView1)
    TNFscalarBar.TitleFontSize = 15
    TNFscalarBar.TitleColor=[0,0,0]
    TNFscalarBar.LabelFontSize=8 
    TNFscalarBar.LabelColor=[0,0,0]
    #TNFscalarBar.DrawSubTickMarks=False
    # TNFscalarBar.MaximumNumberofLabels = 3

    RenderView1.ResetCamera()

    # RenderView1 = GetRenderView()
    # RenderView1.CameraViewUp = [0, 0, 0]
    # RenderView1.CameraPosition = [0, 0, 640]
    # RenderView1.CameraClippingRange = [0, 0]

    camera = GetActiveCamera()
    #camera.Dolly(1.4)

    # camera.SetFocalPoint(0, 0, 0)
    # camera.SetPosition(0, 0, -10)
    # camera.SetViewUp(0, 1, 0)
    # camera.SetViewAngle(30)

    # renderView1 = GetActiveViewOrCreate('RenderView')
    # renderView1.ViewSize = [1670, 1091]

    # outDisplay = GetDisplayProperties(reader, view=renderView1)

    Render()
    
    ExportView(output_TNF_path, view=RenderView1, Plottitle='ParaView GL2PS Export')

    Delete(TNFscalarBar)
    Delete(RenderView1)

    RenderView1 = CreateRenderView()
    RenderView1.UseColorPaletteForBackground=0
    SetViewProperties(Background = [1, 1, 1])
    RenderView1.OrientationAxesVisibility = 0

    Display1=Show()
    camera = GetActiveCamera()

    ##TGBb

    ColorBy(Display1,'TGFb')
    Display1.RescaleTransferFunctionToDataRange(True)

    dataTGFColorMap = GetColorTransferFunction('TGFb')
    dataTGFOpacityMap = GetOpacityTransferFunction('TGFb')

    dataTGFColorMap.EnableOpacityMapping = True

    dataTGFColorMap.ApplyPreset('Oranges', True)
    dataTGFColorMap.InvertTransferFunction()
    
    dataTGFColorMap.RescaleTransferFunction(0.0, 3e-11)
    dataTGFOpacityMap.RescaleTransferFunction(0.0, 3e-11)

    Display1.SetScalarBarVisibility(RenderView1, True)
    TGFscalarBar=GetScalarBar(dataTGFColorMap,RenderView1)
    TGFscalarBar.TitleFontSize = 15
    TGFscalarBar.TitleColor=[0,0,0]
    TGFscalarBar.LabelFontSize=8 
    TGFscalarBar.LabelColor=[0,0,0]
    #TGFscalarBar.DrawSubTickMarks=False
    # TGFscalarBar.MaximumNumberofLabels = 3

    RenderView1.ResetCamera()

    # RenderView1 = GetRenderView()
    # RenderView1.CameraViewUp = [0, 0, 0]
    # RenderView1.CameraPosition = [0, 0, 640]
    # RenderView1.CameraClippingRange = [0, 0]

    camera = GetActiveCamera()
    #camera.Dolly(1.4)

    # camera.SetFocalPoint(0, 0, 0)
    # camera.SetPosition(0, 0, -10)
    # camera.SetViewUp(0, 1, 0)
    # camera.SetViewAngle(30)

    # renderView1 = GetActiveViewOrCreate('RenderView')
    # renderView1.ViewSize = [1670, 1091]

    # outDisplay = GetDisplayProperties(reader, view=renderView1)

    Render()
    
    ExportView(output_TGF_path, view=RenderView1, Plottitle='ParaView GL2PS Export')

    Delete(TGFscalarBar)
    Delete(RenderView1)

    RenderView1 = CreateRenderView()
    RenderView1.UseColorPaletteForBackground=0
    SetViewProperties(Background = [1, 1, 1])
    RenderView1.OrientationAxesVisibility = 0

    Display1=Show()
    camera = GetActiveCamera()

    ##IL10

    ColorBy(Display1,'IL10')
    Display1.RescaleTransferFunctionToDataRange(True)

    # IL10data = GetColorTransferFunction('debris')
    # IL10data.ApplyPreset('Cold and Hot', True)

    dataIL10ColorMap = GetColorTransferFunction('IL10')
    dataIL10OpacityMap = GetOpacityTransferFunction('IL10')

    dataIL10ColorMap.EnableOpacityMapping = True

    dataIL10ColorMap.ApplyPreset('Blues', True)
    dataIL10ColorMap.InvertTransferFunction()
    
    dataIL10ColorMap.RescaleTransferFunction(0.0, 2e-11)
    dataIL10OpacityMap.RescaleTransferFunction(0.0, 2e-11)

    Display1.SetScalarBarVisibility(RenderView1, True)
    IL10scalarBar=GetScalarBar(dataIL10ColorMap,RenderView1)
    IL10scalarBar.TitleFontSize = 15
    IL10scalarBar.TitleColor=[0,0,0]
    IL10scalarBar.LabelFontSize=8 
    IL10scalarBar.LabelColor=[0,0,0]
    #IL10scalarBar.DrawSubTickMarks=False
    # IL10scalarBar.MaximumNumberofLabels = 3

    RenderView1.ResetCamera()

    # RenderView1 = GetRenderView()
    # RenderView1.CameraViewUp = [0, 0, 0]
    # RenderView1.CameraPosition = [0, 0, 640]
    # RenderView1.CameraClippingRange = [0, 0]

    camera = GetActiveCamera()
    #camera.Dolly(1.4)

    # camera.SetFocalPoint(0, 0, 0)
    # camera.SetPosition(0, 0, -10)
    # camera.SetViewUp(0, 1, 0)
    # camera.SetViewAngle(30)

    # renderView1 = GetActiveViewOrCreate('RenderView')
    # renderView1.ViewSize = [1670, 1091]

    # outDisplay = GetDisplayProperties(reader, view=renderView1)

    Render()
    
    ExportView(output_IL10_path, view=RenderView1, Plottitle='ParaView GL2PS Export')

    Delete(IL10scalarBar)
    Delete(RenderView1)

    RenderView1 = CreateRenderView()
    RenderView1.UseColorPaletteForBackground=0
    SetViewProperties(Background = [1, 1, 1])
    RenderView1.OrientationAxesVisibility = 0

    Display1=Show()
    camera = GetActiveCamera()

    ##IFNg

    ColorBy(Display1,'IFNg')
    Display1.RescaleTransferFunctionToDataRange(True)

    # IFNgdata = GetColorTransferFunction('debris')
    # IFNgdata.ApplyPreset('Cold and Hot', True)

    dataIFNgColorMap = GetColorTransferFunction('IFNg')
    dataIFNgOpacityMap = GetOpacityTransferFunction('IFNg')

    dataIFNgColorMap.EnableOpacityMapping = True

    dataIFNgColorMap.ApplyPreset('Purples', True)
    dataIFNgColorMap.InvertTransferFunction()
    
    dataIFNgColorMap.RescaleTransferFunction(0.0, 1e-12)
    dataIFNgOpacityMap.RescaleTransferFunction(0.0, 1e-12)

    Display1.SetScalarBarVisibility(RenderView1, True)
    IFNgscalarBar=GetScalarBar(dataIFNgColorMap,RenderView1)
    IFNgscalarBar.TitleFontSize = 15
    IFNgscalarBar.TitleColor=[0,0,0]
    IFNgscalarBar.LabelFontSize=8 
    IFNgscalarBar.LabelColor=[0,0,0]
    # IFNgscalarBar.MaximumNumberofLabels = 3

    RenderView1.ResetCamera()

    # RenderView1 = GetRenderView()
    # RenderView1.CameraViewUp = [0, 0, 0]
    # RenderView1.CameraPosition = [0, 0, 640]
    # RenderView1.CameraClippingRange = [0, 0]

    camera = GetActiveCamera()
    #camera.Dolly(1.4)

    # camera.SetFocalPoint(0, 0, 0)
    # camera.SetPosition(0, 0, -10)
    # camera.SetViewUp(0, 1, 0)
    # camera.SetViewAngle(30)

    # renderView1 = GetActiveViewOrCreate('RenderView')
    # renderView1.ViewSize = [1670, 1091]

    # outDisplay = GetDisplayProperties(reader, view=renderView1)

    Render()
    
    ExportView(output_IFNg_path, view=RenderView1, Plottitle='ParaView GL2PS Export')

    Delete(IFNgscalarBar)
    Delete(RenderView1)
    #Render()