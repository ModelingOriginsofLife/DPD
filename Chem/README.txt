DPD

Original code:
Thomas Maeke, John McCaskill

Mods for adding chemistry:
Andy Buchanan, Mark Bedau, Norman Packard

Installation:
------------------------
1. Uncompress the dpd package (with, e.g. tar -xzf dpd-????????.tgz)
2. In the DPD directory: make

If you get errors during the make process, check the requirements below.

Requirements:
________________________
DPD requires the SDL 1.2 developer materials: http://www.libsdl.org/
It shouldn't matter how you install them provided that you got the developer libraries and not just the runtime libraries and sdl-config works.

DPD also requires basic software development tools such as gcc, ld and make.  MacOSX and even some modern Linux distributions do not come with these installed by default.  For Mac, downloading the latest XCode software  http://developer.apple.com/tools/xcode/index.html (probably requires free registration) will get you everything you need.  Fink http://fink.sourceforge.net/ is also a good option for getting software.  Linux users should try whatever package-management software came with their distribution if they don't have these.

Usage:
________________________
(from ./dpd3 -h)
   dpd3 [options] [ctrlfiles] [options]

   option:
       -gx        x=geometry
       -tfn       fn=write template ctrl file
       -ffn       fn=ctrl file

       -dx     --disp=x     x=displayinterval [1]
       -h,?    --help       this help
       -v                   verbose
               --dim=x      x=dimension 2 or 3 [3]
               --dt=x.x     x.x=dt, integration time step [0.01]
       -bx.x   --bond=x.x   x.x=bond force between chained monomers [10.0]
       -rx.x   --rho=x.x    x.x=rho,   density of particles per unit cube [4.0]
               --sigma=x.x  x.x=sigma, factor for random and dissipative force [3.0]
               --rfac=x.x   x.x=rfac,  factor for random force [sqrt(3.0)]
       -sx     --size=x     x=system size [16]
               --seed=x     x=random seed [time]
   ESC: Program end

	Many of these parameters can also be changed while DPD is running.  Press 'h' while it's running to see a list of them and which keys effect which parameters.  In general a lower case letter moves the parameter in one direction, while uppercase moves it in the other.
	While these options allow for changing parameters on the fly, DPD is usually run using a control file (*.dat) so that information is better preserved.
