$-----------------------------------------------------------------------------+
$                                                                             |
$   CASE:       M2D_DCI_Implosion_120                                         |
$                                                                             |
$   Version:    211205                                                        |
$                                                                             |
$   History:    from [m2d-cases-7.0.7]/M2D_Laser_110                          |
$                                                                             |
$   Written by: Rafael Ramis and Fuyuan Wu                                    |
$                                                                             |
$   Abstract:   Laser driven DCI implosion                                    |
$                  ALE hydrodynamics with symmetry preservation               |
$                  Radiation transport and heat conduction                    |
$                  Test ray-tracing package                                   |
$                  Filter radiation and laser deposition                      |
$                                                                             |
$-----------------------------------------------------------------------------+





$==============================================================================
$--Main program
entry multi2d()
{
   time                        = 0.25e-9;  $original 15e-9 
   dt_out                      = 50e-12;  
   dt_max                      = 3e-12;    
   dt_init                     = 1e-12;   
$--Create grid definition, initial state and initial time-step
   hydro                       = Capsule();
   hydro.dt                    = dt_init;
   hydro.output                = 'myoutput(1)';

$--construct a structure to store variables for GA optimizations
   struct(hydro,<gaopt>,ga_initialize_fitness());

   mindex    = Arguments#2;
   fileinp   = "inp_":mindex:".dat";
   fileout   = "fit_":mindex:".dat";
   res       = readPulse(fileinp);
   hydro.gaopt.ltimes          = res.times;
   hydro.gaopt.lpowers         = res.powers;
   hydro.gaopt.elaser          = res.elaser;
   if(hydro.gaopt.elaser>80 |hydro.gaopt.elaser<45) {
      ga_output_fitness(fileout,hydro);
      return(0);
   }
   $display(res.times);
   $display(res.powers);
   $display(res.elaser);
   $return(0);

$--Create structure with integration parameters
   control                     = struct();
   struct(control,<hydro>,M20_Hydrodynamics_Control());
   control.hydro.max_cfl       = 0.6;
   control.hydro.max_var       = 0.1;
   control.hydro.kq1           = 1.8;
   control.hydro.kq2           = 0.5;
   control.hydro.qmode         = 5;
   control.hydro.pfactor       = 2;
   control.hydro.spherical     = 1;
   control.hydro.pext_f        = 'pext_f(4)';
   struct(control,<ray>,M20_RayTracing_Control());
   control.ray.rays_f          = 'Laser(2)';
   control.ray.adjust_f        = 'smooth_laser(2)';
   control.ray.wavelength      = 0.0000350;
   struct(control,<transport>,M20_Transport_Control());
   control.transport.Tmaxvar   = 5e2;
   control.transport.reflex    = 1;
   control.transport.nphis     = 2:4:4:2;
   control.transport.freepath  = 'Material_Freepath(4)';
   control.transport.heatcond  = 'Material_Conductivity(4)';
   control.transport.adjust_f  = 'smooth_radia(2)';
   struct(control.transport,<kalpha>,2.5:0.0100);
   struct(control,<plasma>, M20_Plasma_Control());
   struct(control,<project>,M20_Project_Control());
   control.project.var         = 0.6;

$--Perform the simulation
   times                       = (0 ... ceil(time/dt_out))*dt_out;
   hydro1                      = hydro;
   hydroR                      = hydro;
   DEBUG_LINES                 = layers(hydro.p);
   rext                        = hydro.topo.x0[pos(hydro.topo.y0==0)];
   rext                        = rext[sort(rext)];
   while(1){
      control.hydro.pext_a     = [hydro1];
      hydro1                   = minrho(hydro1,1e-5);

      hydro2                   = M20_Step(control,hydro1);

      hydro2                   = ga_calculate_fitness(hydro2);

      $HydroOutputAtTimes(hydro1,hydro2,times);
      if(M20_Project_Test(hydro2,hydroR,control.project)){
         $p("-"[(1 .. 72)>0]);
         $p("   unit  = Project");
         hydro2             = M20_Project(hydro2,hydroR,control.project);
      }
      hydro1                   = hydro2;
      $if(hydro1.dt<1e-15|max(hydro1.TeN)>1e6){
      if(hydro1.dt<1e-15|max(hydro1.TN)>1e6){
         p("-"[(1...72)>0]);
         p("   unit  = abort");
         p(encode("         dt   = %g",hydro1.dt));
         p(encode("         Tmax = %g",max(hydro1.TN)));
         $HydroOutput(hydro2);

         hydro1.gaopt.rhorf    = 0;
         hydro1.gaopt.rhorc    = 0;
         hydro1.gaopt.rho      = 0;
         hydro1.gaopt.ekimp    = 0;
	 ga_output_fitness(fileout,hydro1);

         break;
      }
      hydro1.dt                = min(hydro1.dt,dt_max);
      if(hydro1.time>max(times)){
        ga_output_fitness(fileout,hydro1);
	break;
      }
   }
   p("-"[(1...72)>0]);
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

   return(o);
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

   rhorfmax   = hydro.gaopt.rhorf;
   rhorcmax   = hydro.gaopt.rhorc;
   rhomax     = hydro.gaopt.rho;
   mimp       = hydro.gaopt.mimp;
   ekimp      = hydro.gaopt.ekimp;
   ptmax      = hydro.gaopt.P;
   timax      = hydro.gaopt.T;
   vrmax      = hydro.gaopt.V;

   frac1      = hydro.mass#1/hydro.tmass;
   frac1n     = cell_to_node(topo.tn,topo.ct,topo.np,frac1);
   rhon       = cell_to_node(topo.tn,topo.ct,topo.np,hydro.rho*(frac1>0));
   rhomax     = max(hydro.rho*(frac1>0):rhomax);
   ptmax      = max(hydro.P  *(frac1>0):ptmax);
   timax      = max(hydro.TN*frac1n:timax);
   vrmax      = max(vr*(rhon>10)*frac1n:vrmax);
   nmass1     = hydro.nmass*frac1n*(vx<0);
   mimp       = max(sum(nmass1):mimp);
   ekimp      = max(sum(0.5*nmass1*vr^2):ekimp);

   xn         = hydro.p[ix];
   yn         = hydro.p[iy];
   rn         = hypot(xn,yn);
   rc         = node_to_cell(topo.tn,(topo.ct)*0,rn); 

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
   rhorDT     = sum(hydro.tmass*(rc<16e-4)/rc^2)/(4*3.14159254);
   rhorcmax   = max(rhorDT,rhorcmax);

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

   fprintf(mfid1,"%8.4f %8.4f %8.2f %6.3f %6.2f %6.2f %8.3f %6.2f %6.2f\n", 
           data[0],data[1],data[2],data[3],data[4],data[5],data[6],
	   data[7],data[8]);

   fclose(mfid1);

   DLibera(_file); DLibera(_data);
   return(DCreaNulo());
}
/*C*/
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
$--Set minimum density
local minrho(hydro,minrho)
{
   hydro      = copy(hydro);
   factor     = max(minrho/hydro.rho,1);
   for(i=1;i<=#hydro.mat;i=i+1){
      hydro.mass#i=factor*hydro.mass#i;
   }
   HydroMass(hydro);
   HydroRho(hydro);
   return(hydro);
}
$==============================================================================
$--Extra output 
local myoutput(hydro)
{
   np          = hydro.topo.np;
   iy          = (1 ... np)*2;
   ix          = iy-1;
   vx          = hydro.v[ix];
   vy          = hydro.v[iy];
   vr          = sqrt(vx*vx+vy*vy);
   
                      m2d_write("vr",    vr);
   if(?hydro.ne)      m2d_write("ne",hydro.ne);
   if(?hydro.depo)    m2d_write("delas",hydro.depo);
   if(?hydro.q)       m2d_write("q",hydro.q);
   if(?hydro.lambda)  m2d_write("lambda",hydro.lambda);
   if(?hydro.dexray)  m2d_write("dexray",hydro.dexray);
}
$==============================================================================
$--Find nodes with the same radius
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
$--Smooth laser deposition
smooth_laser(depo,hydro)
{
   depo     = copy(depo);
   lines    = DEBUG_LINES;
   for(i=1;i<=#lines;i=i+1){
      line  = lines#i;
      depo[line] = M20_AzimuthalFilter(depo[line],hydro.nmass[line],
                   dt=0.0001,power=4,cycles=1);
   }
   return(depo);
}
$==============================================================================
$--Smooth radiation and heat flux transport
smooth_radia(depo,hydro)
{
   depo     = copy(depo);
   lines    = DEBUG_LINES;
   for(i=1;i<=#lines;i=i+1){
      line  = lines#i;
      depo[line] = M20_AzimuthalFilter(depo[line],hydro.nmass[line],
                   dt=0.0001,power=4,cycles=1);
   }
   return(depo);
}
$==============================================================================
$--Capsule design based on Temporal et al., Eur. Phys. J. D (2019) 73:5
$--Initial densities modified to have zero pressure at initial time
local Capsule()
{
   $n           = 40:50:10:30;
   n           = 50:100:50:80;
   r           = double(0.012:0.0700:0.0900:0.1300);
   z           = 1:1:1.0:1.01;
   i           = 1:2:3:4;
   ma          = ma_simple_quadrant(n,r,z,i);
   zi          = sum(ma.qid==1|ma.qid==2);
                 ma_divide(ma,1...zi,1:1:1:1);
                 ma_divide(ma,(zi+1)...#ma.qid,1:1:0:0);
   topo        = NewTopo(ma.t,ma.ct,ma.x,ma.y);
   hydro       = Hydro(topo,gmode=2,xsymmetry=1,ntemp=1);
$--Compute center of cells and region flags 'ic', 'iwall', 'iv'
   rnode    = hypot(ma.x,ma.y);
   xn       = ma.y;
   yn       = ma.y;
   xc       = center_of_cells(ma.t,ma.x);
   yc       = center_of_cells(ma.t,ma.y);
   rc       = hypot(xc,yc);

   cone_a  = 50;
   iwall    = ((yc       )> xc*tan(cone_a/45*atan(1)))&
              ((yc-0.0068)<=xc*tan(cone_a/45*atan(1)))&
              (fabs(xc)>=0.0120) &
              (fabs(rc)<=0.0900);
   iDT     = rc>0.0720 & rc<=0.0850 & 
             yc<=rc*sin(cone_a/45*atan(1));	      
   iCH     = rc>0.0850 & rc<=0.0900 & 
             yc<=rc*sin(cone_a/45*atan(1));	      
   iablator = rc>0.0900 & 
             yc<=rc*sin(cone_a/45*atan(1));	      
   iout     = !iwall & !iDT & !iCH & !iablator;

   tempT    = 273/11605;
   rhoDT    = 0.25*iDT+0.008*iout;
   tiDT     = tempT;
   rhoCH    = 1.1*iCH+0.001*iablator;
   tiCH     = tempT*iCH+100*iablator;
   rhoAu    = (19.2)*(iwall);
   tiAu     = tempT;

   matDT    = MaterialDT();
   matCH    = MaterialCH();
   matAu    = MaterialAu();

   HydroAdd(hydro,matDT,rhoDT,tiDT,tiDT);
   HydroAdd(hydro,matCH,rhoCH,tiCH,tiCH);
   HydroAdd(hydro,matAu,rhoAu,tiAu,tiAu);

   HydroEOS(hydro);
   return(hydro);
}
$==============================================================================
$ set laser pulse shape
local Laser(hydro,dt)
{
   time        = hydro.time;

   $R5-V1
   $times       =     0:  0.85:  1.50:   2.1:  2.5:   2.7: 2.9: 4.753: 4.754;
   $powers      = 0.305: 0.398: 0.585: 1.246: 3.24: 6.919:15.4:  15.4: 0; 

   $R6-V1
   $times  = 0: 0.20: 0.6: 0.8: 1.19: 1.60: 1.90: 2.24: 2.67: 3.12:  3.39: 4.37: 4.67;
   $powers = 0: 1.23: 0.0: 0.0: 1.53: 1.39: 3.47: 5.90: 5.98: 8.63: 11.27: 15.4: 0; 

   $R7-V3
   $times  = 0: 0.13:0.26:0.76: 0.89: 1.02: 1.25: 1.38: 1.51: 1.73:  1.86: 4.81: 4.94;
   $powers = 0: 1.81: 0.0: 0.0: 3.31: 0.00: 0.00: 6.08: 0.00: 0.00:  12.4: 12.4: 0; 

   $120kJDT
   times  = 0: 0.88:1.88:2.68: 3.35: 4.67: 5.42: 5.98: 6.64: 7.08: 7.36 :10.36: 10.49;
   powers = 0: 0.72: 0.0:0.72: 0.34: 0.15: 1.45: 1.87: 3.59: 9.03: 14.29:14.29: 0; 

   times       = hydro.gaopt.ltimes;
   powers      = hydro.gaopt.lpowers;
   power       = algorithm_interpola_lin1d(times,powers,time*1e9,0)*1e19;
   radius      = 0.80*0.0900*sin(50/45*atan(1));
   n           = 200;
   zero        = (1 ... n)*0;
   angle       = 50*atan(1)/45;
   r           = radius*((((0 ... (n-1))+random(n))/n)*2-1);
   a           = random(n)*atan(1)*8;
   x0          = (10*cos(angle)+r*cos(a)*sin(angle))+0.8*0.0900;
   y0          = (10*sin(angle)-r*cos(a)*cos(angle));
   z0          = (r*sin(a));
   x1          = (zero-1*cos(angle));
   y1          = (zero-1*sin(angle));
   z1          = (zero);
   w           = (fabs(r)/exp((fabs(r)/radius)^4));

   w           = power*w/sum(w);
   return([x0,x1,y0,y1,z0,z1,w]);
}
$==============================================================================
$==============================================================================
$--Compute radial distribution of cells
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
$--auxiliar routine
local ph(time,h)
{
   save=h.time;
   h.time=time;
   HydroOutput(h);
   h.time=save;
}
$==============================================================================
local center_of_cells(tn,x)
{
   nt          = #tn/3;
   i1          = tn[3*(2*(1 ... nt/2)-1)-2];
   i2          = tn[3*(2*(1 ... nt/2)-1)-1];
   i3          = tn[3*(2*(1 ... nt/2)-1)  ];
   i4          = tn[3*(2*(1 ... nt/2)-1)+1];
   xc          = (1 ... nt)*double(0);
   xc[2*(1 ... nt/2)-1] = (x[i1]+x[i2]+x[i3]+x[i4])/4.0;
   xc[2*(1 ... nt/2)  ] = (x[i1]+x[i2]+x[i3]+x[i4])/4.0;
   return(xc);
}
$==============================================================================
$--External pressure on the boundary wall
local pext_f(time,xm,ym,hydro)
{  
   bw          = hydro.topo.bw;
   wt          = hydro.topo.wt;
   index2      = 2*bw;
   index1      = index2-1;
   t1          = wt[index1];
   t2          = wt[index2];
   pres        = 0.5*hydro.P[t1+t2];
   return(pres);
}
