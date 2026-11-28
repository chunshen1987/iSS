c     C-interoperable bridge to the UrQMD 3.4 COMMON blocks.
      module iss_urqmd_bridge
      use, intrinsic :: iso_c_binding, only: c_int,c_double
      implicit none

      type, bind(c) :: particle
        real(c_double) :: t,x,y,z,energy,px,py,pz,mass
        real(c_double) :: decay_time,formation_time,cross_section_factor
        integer(c_int) :: ityp,iso3,charge
      end type

      contains

      subroutine iss_urqmd_setup() bind(c)
      include 'coms.f'
      include 'options.f'

c     Initialize UrQMD once, with default options and no burner file.
      eos=0
      call uinit(1)
c     Match the supplied uqmd.burner setting for total cross sections.
      CTParam(22)=1.d0
c     Collision and decay routines still call the UrQMD output entries.
      bf13=.true.
      bf14=.true.
      bf15=.true.
      bf16=.true.
      bf17=.true.
      bf18=.true.
      bf19=.true.
      bf20=.true.
      bf30=.true.
      end subroutine

      integer(c_int) function iss_urqmd_run(count,initial,seed,
     &     steps,step_size,index,collisions) bind(c)
      include 'coms.f'
      include 'comres.f'
      include 'colltab.f'
      include 'newpart.f'
      integer(c_int),value :: count,seed,steps,index
      real(c_double),value :: step_size
      type(particle),intent(in) :: initial(*)
      integer(c_int),intent(out) :: collisions
      integer i

      iss_urqmd_run=-1
      if(count.le.0.or.count.gt.nmax) return
      call iss_urqmd_reset_event
      nct=0
      nsav=0
      actcol=0
      nspec=0
      nexit=0
      event=index
      ranseed=seed
      call sseed(ranseed)
      nsteps=steps
      dtimestep=step_size
      npart=count

      do i=1,npart
        r0(i)=initial(i)%t
        rx(i)=initial(i)%x
        ry(i)=initial(i)%y
        rz(i)=initial(i)%z
        p0(i)=initial(i)%energy
        px(i)=initial(i)%px
        py(i)=initial(i)%py
        pz(i)=initial(i)%pz
        fmass(i)=initial(i)%mass
        dectime(i)=initial(i)%decay_time
        tform(i)=initial(i)%formation_time
        xtotfac(i)=initial(i)%cross_section_factor
        ityp(i)=initial(i)%ityp
        iso3(i)=initial(i)%iso3
        charge(i)=initial(i)%charge
        origin(i)=99
        uid(i)=i
        if(abs(ityp(i)).le.maxbar) nbar=nbar+1
      enddo
      uid_cnt=npart
      nmes=npart-nbar
      acttime=r0(1)
      time=acttime
      call iss_urqmd_cascade(collisions)
      iss_urqmd_run=npart
      end function

      subroutine iss_urqmd_particles(final) bind(c)
      include 'coms.f'
      type(particle),intent(out) :: final(*)
      integer i
      do i=1,npart
c       Export the final phase-space point, as in the first nine
c       fields of file 13 (including the Fermi momentum correction).
        final(i)=particle(r0(i),rx(i),ry(i),rz(i),p0(i),
     &      px(i)+ffermpx(i),py(i)+ffermpy(i),pz(i)+ffermpz(i),
     &      fmass(i),dectime(i),tform(i),xtotfac(i),
     &      ityp(i),iso3(i),charge(i))
      enddo
      end subroutine
      end module
