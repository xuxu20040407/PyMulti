$-----------------------------------------------------------------------------+
$                                                                             |
$   CASE:       M2D_Laser_Fusion                                              |
$                                                                             |
$   Version:    200803                                                        |
$                                                                             |
$   History:    from [m2d-cases-7.0.8]/M2D_Laser_ALE_120                      |
$                                                                             |
$   Written by: Rafael Ramis and Zheyi Ge                                     |
$                                                                             |
$   Abstract:   Laser driven capsule implosion                                |
$                  ALE hydrodynamics with symmetry preservation               |
$                  Radiation transport and heat conduction                    |
$                  Test ray-tracing package                                   |
$                  Filter radiation and laser deposition                      |
$                  DT nuclear reations and alpha-particle transport           |
$                                                                             |
$   Note:       Capsule design based on Temporal et al., Eur. Phys. J. D      |
$               (2019) 73:5. Value obtainde by this simulation                |
$               (absorbed laser 400 kJ and nuclear yield 10 MJ) are similar   |
$               to the 1D published values (270 kJ and 11.9 MJ)               |
$                                                                             |
$-----------------------------------------------------------------------------+





$==============================================================================
$--main program
entry multi2d()
{
$--capsule design based on Temporal et al., Eur. Phys. J. D (2019) 73:5
$--initial densities modified to have zero pressure at initial time
   $r           = double(0.0593:0.0791:0.0815);  $--radius of interfaces
   r           = double(0.0200:0.0450:0.0550);  $--radius of interfaces
   n           = 16:30:20;                      $--radial number of cells
   z           = 1:1.0:0.9;                     $--zone parameters
   i           = 1:2:3;                         $--zone indices
$--use special routine in [malla-6.0] to generate an spherical capsule
   ma          = ma_simple_quadrant(n,r,z,i);
$--use different quadrangle subdivision mode in each layer
   zi          = sum(ma.qid==1|ma.qid==2);       
                 ma_divide(ma,1...zi,1:1:1:1);
                 ma_divide(ma,(zi+1)...#ma.qid,1:1:0:0);
$--generate the main data structure
   topo        = NewTopo(ma.t,ma.ct,ma.x,ma.y);
   h           = Hydro(topo,gmode=2,xsymmetry=1,ntemp=2);
$--three materials DT, He (reaction product), and CH (ablator)
   temp        = 0.01;
   rhoGas      = 0.0005; 
   rhoIce      = 0.0005; 
   rhoCH       = 1.127320;

   $matDT       = MaterialDT();
   $matHe       = MaterialDT();

   matCH       = MaterialCH();
   HydroAdd(h,matCH,(rhoGas:rhoIce:rhoCH)[ma.tid],temp,temp);

$--groups of nodes with identical initial radius
   lines       = layers(h.p);
$--create structure with integration parameters
   control                     = struct();
   struct(control,<hydro>,    M20_Hydrodynamics_Control());
   struct(control,<ray>,      M20_RayTracing_Control());
   struct(control,<project>,  M20_Project_Control());
   struct(control,<transport>,M20_Transport_Control());
   struct(control,<plasma>,   [[]]);
   struct(control,<fusion>,   [[<f>]]);
   control.hydro.max_cfl       = 1;
   control.hydro.max_var       = 0.1;
   control.hydro.qmode         = 5;
   control.hydro.pfactor       = 2;
   control.hydro.spherical     = 1;
   control.ray.rays_f          = 'Laser(2)';
   control.ray.adjust_f        = 'smooth_laser(3)';
   control.ray.adjust_a        = [lines];
   control.ray.wavelength      = 0.0000350;
   control.transport.Tmaxvar   = 1e6;
   control.transport.reflex    = 1;
   control.transport.nphis     = 2:4:4:4:4:2;
   control.transport.freepath  = 'Material_Freepath(4)';
   control.transport.heatcond  = 'Material_Conductivity(4)';
   control.transport.adjust_f  = 'smooth_radia(3)';
   control.transport.adjust_a  = [lines];
   struct(control.transport,<kalpha>,1.5:0.0010);
$--additional output    
   struct(h,<output>,'myoutput(1)');
$--initial time step
   h.dt        = 10e-12;
$--output at these times
   times       =                 (0 ... 100)*0.100e-9;
$--initial radius of nodes on the horizontal axis
   rext        = h.topo.x0[pos(h.topo.y0==0)];
   rext        = rext[sort(rext)];


$--construct a structure to store variables for GA optimizations
   struct(h,<gaopt>,ga_initialize_fitness());
   
   mindex    = Arguments#2;
   fileinp   = "inp_":mindex:".dat";
   fileout   = "fit_":mindex:".dat";
   res       = readPulse(fileinp);
   h.gaopt.ltimes          = res.times;
   h.gaopt.lpowers         = res.powers;
   h.gaopt.elaser          = res.elaser;
   if(h.gaopt.elaser>1000 |h.gaopt.elaser<0) {
      h.gaopt.amp = 0.01;
      ga_output_fitness(fileout,hydro);
      return(0);
   }

$--values at the beginning of time step
   h1          = h;   
$--reference grid
   hR          = h;   

   while(1){
$-----perform numerical integration
      h2       = M20_Step(control,h1);

      h2       = ga_calculate_fitness(h2);

$-----print interesting values 
      rhor(h2);
$-----output at specified times
      $HydroOutputAtTimes(h1,h2,times);

$-----after 5 ns, consider a grid change
      if(M20_Project_Test(h2,hR,control.project)&h2.time>5e-9){
         p("-"[(1 .. 72)>0]);
         p("   unit  = Project");
$--------try to obtain a smooth grid based of interface #47 and time
         hb    = procesa(h2,rext[57],h2.time>6.0e-9);
         $hb    = procesa(h2,rext[47],h2.time>6.0e-9);
$--------if successful, project values and take the new grid as reference
         if(?hb){
            p("   perform projection");
            h2 = M20_Project(h2,hb,control.project);
            hR = h2;
         }
      }
$-----check that time step and maximum temperature are inside bounds
      if(h2.dt<1e-15|max(h2.TeN)>1e6){
         p("-"[(1...72)>0]);
         p("   unit  = abort");
         p(encode("         dt   = %g",h2.dt));
         p(encode("         Tmax = %g",max(h2.TeN)));
         HydroOutput(h2);

         ga_output_fitness(fileout,h2);

         break;
      }
$-----normal exit condition 
      if(h2.time>max(times)){
         ga_output_fitness(fileout,h2);

         break;
      }

$-----step finished, replace initial values for next step
      h1       = h2;
$-----set maximum time step
      h1.dt    = min(h1.dt,10e-12);
   }
   p("-"[(1...72)>0]);
}
$==============================================================================
$--specific output for thermonuclear reactions
local myoutput(hydro){

   np          = hydro.topo.np;
   iy          = (1 ... np)*2;
   ix          = iy-1;
   vx          = hydro.v[ix];
   vy          = hydro.v[iy];
   topo        = hydro.topo;
   rhon        = cell_to_node(topo.tn,topo.ct,topo.np,hydro.rho);
   vr          = sqrt(vx*vx+vy*vy)*(vx<0)*(rhon>0.1);
   
                      m2d_write("vr",    vr);
   if(?hydro.mfp_alpha)   m2d_write("mfp_alpha",hydro.mfp_alpha);
   if(?hydro.fusion)      m2d_write("fusion",   hydro.fusion);
   if(?hydro.alpha_dep)   m2d_write("alpha_dep",hydro.alpha_dep);
}
$==============================================================================
$--compute and print confinement parameters
local rhor(h)
{
   a           = algorithm_components(2,h.p);
   r           = hypot(a#1,a#2);
   a           = algorithm_components(3,h.topo.tn);
   rm          = (r[a#1]+r[a#2]+r[a#3])/3;
   rhor        = sum(((h.mass)#1)/rm^2)/(4*3.14159254);
   p("   confinement parameters");
   p(encode("   rhor   = %g",rhor));
   rhorHS      = sum((h.Te>5000)*((h.mass)#1)/rm^2)/(4*3.14159254);
   p(encode("   rhorHS = %g",rhorHS));
}
$==============================================================================
$--find groups of nodes with the same radius
local layers(p)
{
   a           = algorithm_components(2,p);
   x           = a#1;
   y           = a#2;
   r           = hypot(x,y);
   i           = sort(r);
   r           = r[i];
   dr          = cut_first(r)-cut_last(r);
   j1          = pos(1:(dr>max(r)*1e-6));
   j2          = (cut_first(j1)-1):#r;
   ind         = pos(j2-j1==max(j2-j1));
   j1          = j1[ind];
   j2          = j2[ind];
   lines       = split(i,j1,j2); 
   for(i=1;i<=#lines;i=i+1){
      line     = lines#i;
      ind      = sort(atan2(y[line],x[line]));
      lines#i  = line[ind];
   }
   return(lines);
}
$==============================================================================
$--smooth laser deposition over groups of nodes with same initial radius
smooth_laser(depo,hydro,lines)
{
   depo     = copy(depo);
   for(i=1;i<=#lines;i=i+1){
      line  = lines#i;
      depo[line] = M20_AzimuthalFilter(depo[line],hydro.nmass[line],
                   dt=0.0001,power=4,cycles=1);
   }
   return(depo);
}
$==============================================================================
$--smooth radiation and heat flux over groups of nodes with same initial radius
smooth_radia(depo,hydro,lines)
{
   depo     = copy(depo);
   for(i=1;i<=#lines;i=i+1){
      line  = lines#i;
      depo[line] = M20_AzimuthalFilter(depo[line],hydro.nmass[line],
                   dt=0.0002,power=4,cycles=1);
   }
   return(depo);
}
$==============================================================================
local Laser(hydro,dt)
{
   time        = hydro.time;

   times       = hydro.gaopt.ltimes;
   powers      = hydro.gaopt.lpowers;   

   power       = algorithm_interpola_lin1d(times,powers,time*1e9,0)*1e19;
   power       = power*0.5;

   n           = 160;
   zero        = (1 .. n)*0;
   one         = zero+1;
   theta       = ((1 .. n)-1+0.5)/(n)*atan(1)*2;
   x0          = 2*cos(theta);
   y0          = 2*sin(theta);
   z0          = 1e-6;
   x1          = -1*cos(theta);
   y1          = -1*sin(theta);
   z1          = -0.5e-6;
$--Supergaussian profile
   w          = y0;
$--Rescale beamlet power
   w           = power*w/sum(w);
   return([x0,x1,y0,y1,z0,z1,w]);

}
$==============================================================================
$--regriding package (currently under development)
local procesa(h,rext,flag)
{
   h           = copy(h);
   h.p         = copy(h.p);
$--Number pi     
   pi          = 4*atan(double(1));
$--Small quantity
   eps         = 1e-6;
$--Obtain region to be modified
   x0          = h.topo.x0;
   y0          = h.topo.y0;
   r0          = hypot(x0,y0);
   nodes       = pos(r0>0&r0<rext*(1-eps));
   a0          = atan2(y0[nodes],x0[nodes])/(pi/2);
   a0[pos(fabs(a0-1)<eps)] = 1;
$--Obtain exterior profile
   extnodes    = pos(r0<rext*(1+eps)&r0>rext*(1-eps));
   extnodes    = extnodes[sort(atan2(y0[extnodes],x0[extnodes]))];
$--Circularize exterior profile
   if(flag){
      xcir     = h.p[extnodes+extnodes-1];
      ycir     = h.p[extnodes+extnodes];
      rcir     = hypot(xcir,ycir);
      acir     = atan2(ycir,xcir);
      acir     = (0 ... (#acir-1))/(#acir-1)*max(acir);
      rcir     = sum(rcir)/#rcir;
      xcir     = rcir*cos(acir);
      ycir     = rcir*sin(acir);
      xcir[#xcir]=0;
      h.p[extnodes+extnodes-1]=xcir;
      h.p[extnodes+extnodes]=ycir;
   }
$--Compute external n radius and azimuth
   x0b         = x0[extnodes];
   y0b         = y0[extnodes];
   a0b         = atan2(y0b,x0b)/(pi/2);
   a0b[pos(fabs(a0b-1)<eps)] = 1;
   x1b         = h.p[extnodes+extnodes-1];
   y1b         = h.p[extnodes+extnodes];
   r1b         = hypot(x1b,y1b);
   a1b         = atan2(y1b,x1b)/(pi/2);
   a1b[pos(fabs(a1b-1)<eps)] = 1;
$--Report excesive distortion
   da          = cut_first(a1b)-cut_last(a1b);
   dist        = min(da)/(sum(da)/#da);
   if(dist<0.5){
      p("dist<0.5");
   }
   db          = (max(r1b)-min(r1b))/max(r1b);
   if(db>0.25){
      p("db<0.25");
   }
$--Compute radial distribution of cells
   map         = map(h,rext); 
$--Interpolate radius and azimuth 
   r0          = r0[nodes];
   factor      = algorithm_interpola_lin1d(map.r0,map.f,r0,r0);
   r1          = algorithm_interpola_lin1d(a0b,r1b,a0,a0)*factor;
   a1          = algorithm_interpola_lin1d(a0b,a1b,a0,a0)*(pi/2);
   s1          = sin(a1); 
   c1          = cos(a1);
   c1[pos(fabs(c1)<eps)] = 0;
   x1          = r1*c1;
   y1          = r1*s1;
$--Return new grid
   h.p[nodes+nodes-1] = x1;
   h.p[nodes+nodes]   = y1;
   HydroVol(h);
   if(min(h.area)<0){
      HydroOutput(h);
      p("min(h.area)<0");
      exit(0);
   }
   return(h);
}
$==============================================================================
$--compute radial distribution of cells (currently under development)
local map(h,rext)
{
   eps         = 1e-6;
   x0          = h.topo.x0;
   y0          = h.topo.y0;
   r0          = hypot(x0,y0);
   r0c         = node_to_cell(h.topo.tn,h.topo.ct,r0);
   latnodes    = pos(r0<rext*(1+eps)&y0==0);
   latnodes    = latnodes[sort(r0[latnodes])];
   r0lat       = h.topo.x0[latnodes];
   r1lat       = h.p[latnodes+latnodes-1];
   vlat        = double(latnodes*0);
   mlat        = double(latnodes*0);
   for(i=1;i<#latnodes;i=i+1){
      a1       = r0lat[i];
      a2       = r0lat[i+1];
      sel      = pos(r0c>=a1&r0c<=a2);
      f1       = (a2-r0c[sel])/(a2-a1);
      f2       = 1-f1;
      vol      = h.vol[sel];
      mass     = vol*(h.rho[sel]);
      mlat[i]  = mlat[i]+sum(f1*mass);
      mlat[i+1]= mlat[i+1]+sum(f2*mass);
      vlat[i]  = vlat[i]+sum(f1*vol);
      vlat[i+1]= vlat[i+1]+sum(f2*vol);
   }
   rholat      = mlat/vlat;
   rholat      = sqrt(rholat);
   ilat        = 0:acc((cut_first(rholat)+cut_last(rholat))*
                       (cut_first(r1lat)-cut_last(r1lat)));
   flat        = (r1lat-min(r1lat))/(max(r1lat)-min(r1lat));
   ilat        = ilat/max(ilat);
   ilat        = ilat/max(ilat)+0.5*flat;
   ilat        = ilat/max(ilat);
   i0lat       = 0 ... (#latnodes-1);
   i0lat       = i0lat/max(i0lat);
   f           = algorithm_interpola_lin1d(ilat,flat,i0lat,i0lat);
   map         = struct();
   struct(map,<r0>,r0lat);
   struct(map,<f>,f);
   return(map);
}
$==============================================================================
local ga_initialize_fitness()
{
   o          = struct();
   struct(o,<rhorf>,   0);
   struct(o,<rhorc>,   0);
   struct(o,<rho>,     0);
   struct(o,<T>,       0);
   struct(o,<P>,       0);
   struct(o,<V>,       0);
   struct(o,<ltimes>,  *);
   struct(o,<lpowers>, *);
   struct(o,<mimp>,    0);
   struct(o,<ekimp>,   0);
   struct(o,<elaser>,  0);
   struct(o,<amp>,     0);
   return(o);
}
$==============================================================================
$ Read input laser pulse  
local readPulse(file)
{
   list     = separa(carga(file),"&");
   $all      = double(decode(list#1));
   all      = decode(list#1);
   n        = (#all-3)/2;
   dt       = all[  1 ...   n];    
   times    = acc(dt);
   powers   = all[n+1 ... 2*n];

   elaser   = 0.5*sum((powers[1 ...n-1]+powers[2 ...n])*dt[2 ...n]);   

   res      = [[<times>,<powers>,<elaser>],times,powers,elaser];
   return(res);
}
$==============================================================================
local ga_calculate_fitness(hydro)
{
   hydro      = copy(hydro);

   topo       = hydro.topo;

   np         = hydro.topo.np;
   iy         = (1 ... np)*2;
   ix         = iy-1;
   vx         = hydro.v[ix];
   vy         = hydro.v[iy];
   vr         = sqrt(vx*vx+vy*vy)*(vx<0);

   xn         = hydro.p[ix];
   yn         = hydro.p[iy];
   rn         = hypot(xn,yn);
   rc         = node_to_cell(topo.tn,(topo.ct)*0,rn); 

   rhorfmax   = hydro.gaopt.rhorf;
   rhorcmax   = hydro.gaopt.rhorc;
   rhomax     = hydro.gaopt.rho;
   mimp       = hydro.gaopt.mimp;
   ekimp      = hydro.gaopt.ekimp;
   ptmax      = hydro.gaopt.P;
   timax      = hydro.gaopt.T;
   vrmax      = hydro.gaopt.V;

   frac2      = hydro.mass#1/hydro.tmass;
   frac2n     = cell_to_node(topo.tn,topo.ct,topo.np,frac2);
   rhon       = cell_to_node(topo.tn,topo.ct,topo.np,hydro.rho*(frac2>0));
   rhomax     = max(hydro.rho*(frac2>0):rhomax);
   ptmax      = max(hydro.Pe  *(frac2>0):ptmax);
   timax      = max(hydro.TeN*(rn>0.0050):timax);
   vrmax      = max(vr*(rhon>1e-1):vrmax);
   nmass1     = hydro.nmass*frac2n*(vx<0);
   mimp       = max(sum(nmass1):mimp);
   ekimp      = max(sum(0.5*nmass1*vr^2):ekimp);

   ind_yaxis  = pos(xn ==0);
   rn_yaxis   = yn[ind_yaxis];
   rhon_yaxis = rhon[ind_yaxis];
   ind_yaxis  = sort(rn_yaxis);
   rn_yaxis   = rn_yaxis[ind_yaxis];
   rhon_yaxis = rhon_yaxis[ind_yaxis];
   dr_yaxis   = cut_first(rn_yaxis)-cut_last(rn_yaxis);
   rhoc_yaxis = 0.5*(cut_first(rhon_yaxis)+cut_last(rhon_yaxis));
   rhorfmax   = max(sum(rhoc_yaxis*dr_yaxis),rhorfmax);
   
   ind_xaxis  = pos(yn == 0);
   rn_xaxis   = xn[ind_xaxis];
   rhon_xaxis = rhon[ind_xaxis];
   ind_xaxis  = sort(rn_xaxis);
   rn_xaxis   = rn_xaxis[ind_xaxis];
   rhon_xaxis = rhon_xaxis[ind_xaxis];
   dr_xaxis   = cut_first(rn_xaxis)-cut_last(rn_xaxis);
   rhoc_xaxis = 0.5*(cut_first(rhon_xaxis)+cut_last(rhon_xaxis));
   rhorcmax   = max(sum(rhoc_xaxis*dr_xaxis),rhorcmax);

   $rhorDT     = sum(((hydro.mass)#1)*(rc<16e-4)/rc^2)/(4*3.14159254);
   $rhorDT     = sum(hydro.tmass*(rc<16e-4)/rc^2)/(4*3.14159254);
   $rhorcmax   = max(rhorDT,rhorcmax);

   hydro.gaopt.rhorf = rhorfmax;
   hydro.gaopt.rhorc = rhorcmax;
   hydro.gaopt.rho   = rhomax;
   hydro.gaopt.P     = ptmax;
   hydro.gaopt.T     = timax;
   hydro.gaopt.V     = vrmax;
   hydro.gaopt.mimp  = mimp;
   hydro.gaopt.ekimp = ekimp;
   return(hydro);
}

$==============================================================================

local ga_output_fitness(fileout,hydro)
{ 
  rhorf      = hydro.gaopt.rhorf;         $ g/cm^2
  rhorc      = hydro.gaopt.rhorc;         $ g/cm^2
  rho        = hydro.gaopt.rho;           $ g/cc
  vr         = hydro.gaopt.V*1e-5;        $ km/s
  ti         = hydro.gaopt.T*1e-3;        $ keV
  ptmax      = hydro.gaopt.P*1e-15;       $ erg/cc-->Gbar
  mimp       = hydro.gaopt.mimp*1e6;      $ ug
  ekimp      = hydro.gaopt.ekimp*1e-10;   $ kJ
  elaser     = hydro.gaopt.elaser;        $ kJ
  amp        = hydro.gaopt.amp;           $ cm

  $dataout    = elaser:amp;
  dataout    = rhorf:rhorc:rho: ti:vr:ptmax: mimp:ekimp:elaser;
  ga_writeVector(fileout,double(dataout));
}
$==============================================================================
$ data -A vector to be written with a type of  double
local ga_writeVector(file,data);

/*C*/
#include <stdio.h>
#include <string.h>
static D *_ga_writeVector2(D* _file, D* _data)
{
   char   s[80]; 
   double *data;
   FILE   *mfid1;

   data = _data->p.d;

   sprintf(s,"%s",_file->p.c);
   mfid1  = fopen(s,"w");

   fprintf(mfid1,"%8.5f %8.5f %8.2f  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.2f\n", 
           data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8]);

   fclose(mfid1);

   DLibera(_file); DLibera(_data);
   return(DCreaNulo());
}
/*C*/
$==============================================================================
