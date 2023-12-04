# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

sim = input('Specify simulation folder')
x_d = float(input('Specify the start position'))
x_pos = 0.33533+(x_d*0.0064)
tube_length = 0.91
end = x_pos + tube_length
y_pos = float(input('specify radial midpoint coords for y-centroid'))
z_pos = float(input('specify centripidal midpoint coords for z-centroid'))
end_time = float(input('Specify end time'))

# get active source.
pIVfoam = GetActiveSource()

# set active source
SetActiveSource(pIVfoam)

UpdatePipeline(time=end_time, proxy=pIVfoam)

# get animation scene
animationScene1 = GetAnimationScene()

# Properties modified on animationScene1
animationScene1.AnimationTime = end_time

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=pIVfoam)
plotOverLine1.Point1 = [x_pos, y_pos, z_pos]
plotOverLine1.Point2 = [end, y_pos, z_pos]

UpdatePipeline(time=end_time, proxy=plotOverLine1)

# save data
SaveData('C:/Users/drewm/Documents/2.0 MSc/2.0 Simulations/' + sim + '/axial_data.csv', proxy=plotOverLine1, ChooseArraysToWrite=1,
    PointDataArrays=['Ma', 'MaT', 'kappa', 'R', 'RMean', 'T', 'TMean', 'TPrime2Mean', 'U', 'UMean', 'UPrime2Mean', 'alphat', 'epsilon', 'k', 'kMean', 'kPrime2Mean', 'nut', 'nutMean', 'nutPrime2Mean', 'p', 'pMean', 'pPrime2Mean', 'rho', 'rhoMean', 'rhoPrime2Mean', 'turbulenceProperties:R', 'turbulenceProperties:omega',])

# set active source
SetActiveSource(pIVfoam)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOverLine1)

# destroy plotOverLine1
Delete(plotOverLine1)
del plotOverLine1