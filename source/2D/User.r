$-----------------------------------------------------------------------------+
$                                                                             |
$   CASE:       M2D_Implo                                                     |
$                                                                             |
$   Version:    200513                                                        |
$                                                                             |
$   History:    from [m2d-cases-7.0.6]/M2D_Implo_100                          |
$                                                                             |
$   Written by: Rafael Ramis                                                  |
$                                                                             |
$   Abstract:   Check symmetry preserving in Lagrangian hydrodynamics         |
$                                                                             |
$-----------------------------------------------------------------------------+





$==============================================================================
$--Routine with user defined data
local problem()
{
   problem=[[<r0>,<dr>,<rho_core>,<rho_dt>,<amp>,<num>,<legendre2>,
   <legendre4>,<res>]];
$--Geometry
   problem.r0        = 0.2;
   problem.rho_core  = 0.0003;
   problem.dr        = 0.00864;
   problem.rho_dt    = 0.220533;
$--Relative resolution
   problem.res       = 1;
$--Initial distortion r=R*(1+amp*cos(theta*num)
   problem.amp       = 0.001;
   problem.num       = 12;         
$--Pressure nonuniformity
   problem.legendre2 = 0;
   problem.legendre4 = 0;
   return(problem);
}
$==============================================================================
$--Main program
entry multi()
{
$--Get users setable parameters
   problem     = problem();
$--Generate grid 
   ma          = ma();
   n1          = round(20*problem.res); 
   n2          = round(30*problem.res);
   nr          = round(20*problem.res);
   r1          = problem.r0;
   r2          = problem.r0+problem.dr;
   p0          = ma_point(ma,0,0);
   p1          = ma_point(ma,r1,0);
   p2          = ma_point(ma,0,r1);
   p3          = ma_point(ma,r2,0);
   p4          = ma_point(ma,0,r2);
   p5          = ma_point(ma,2*r1/n1,0);
   p6          = ma_point(ma,0,2*r1/n1);
                 ma_line(ma, p5, p6,   ma_d(nr),2*r1/n1);
                 ma_line(ma, p5, p1,   ma_d(n1-1,1),0);
                 ma_line(ma, p1, p2,   ma_d(nr),r1);
                 ma_line(ma, p1, p3,   ma_d(n2,1),0);
                 ma_line(ma, p3, p4,   ma_d(nr),r2);
$--Fill regions (1-gas, 2-solid)
   q           = ma_rect(ma,p1,p2,p6,p5,1); 
                 ma_divide(ma,q,0:1:1:0);
   q           = ma_rect(ma,p1,p2,p4,p3,2); 
                 ma_divide(ma,q,1:0:0:1);
                 ma_fill_cap(ma,p0,ma_getline(ma,p5,p6),1);
$--Distort the target
   x           = ma.x;
   y           = ma.y;
   f           = 1+problem.amp*cos(problem.num*atan2(ma.y,ma.x));
   ma.x        = ma.x*f;
   ma.y        = ma.y*f;
$--Generate topology and hydro structure
   m           = ma.tid;
   topo        = NewTopo(ma.t,ma.ct,ma.x,ma.y);
   hydro       = Hydro(topo,gmode=2,xsymmetry=1,temp=1);
   hydro.dt    = 1e-12;
   struct(hydro,<output>,'myoutput(1)');
$--Obtain material (DT with tabulated equation-of-state)
   matDT        = MaterialDT();
   HydroAdd(hydro,matDT,(problem.rho_core:problem.rho_dt)[ma.tid],temp=0,*);
$--Create simulation state and control structures
   control                     = struct();
   struct(control,<hydro>,M20_Hydrodynamics_Control());
   control.hydro.max_cfl       = 1;
   control.hydro.max_var       = 0.5;
   control.hydro.kq1           = 1;
   control.hydro.pext_f        = 'external(4)';
   control.hydro.pext_a        = [problem];
   control.hydro.spherical     = 1;
   control.hydro.pfactor       = 1.2;
   control.hydro.qmode         = 5;
$--Simulation loop
   hydro1                      = hydro;
   times                       = (0 ... 120)*1e-10;
   while(1){
      hydro2=M20_Step(control,hydro1);
      HydroOutputAtTimes(hydro1,hydro2,times);
      if(hydro2.time>=max(times))break;
      if(hydro2.dt<1e-21){
         p("-"[(1 .. 72)>0]);
         p("   dt too small");
         break;
      }
      hydro1=hydro2;
   }
   p("-"[(1 .. 72)>0]);
}
$==============================================================================
$--Additional output
local myoutput(hydro)
{
   if(?hydro.cfl)m2d_write("cfl",hydro.cfl);
   if(?hydro.var)m2d_write("var",hydro.var);
}
$==============================================================================
$--External pressure
local external(time,xs,ys,problem)
{
   if     (time<2.6e-9)   pmean=1.e12;
   else if(time<3.2e-9)   pmean=4.e12;
   else                   pmean=16.e12;
   mu=xs/hypot(xs,ys);
   l2=(3*mu^2-1)/2;
   l4=(35*mu^4-30*mu^2+3)/8;
   f=1+problem.legendre2*l2+problem.legendre4*l4;
   return(pmean*f);
}
$==============================================================================
