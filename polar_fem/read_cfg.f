c ====== Read configuration file
      subroutine read_cfg()
        include 'th.fh'
c === set default values
        DATA DIM_L/1.0/, DIM_D/1.0/, DIM_Lw/0.3/, DIM_Dw/0.5/,
     &       DIM_Ld/0.01/, DIM_Dd/0.02/, DIM_Dr/0.0/,
     &       DIM_Li/1.0/, DIM_Lo/1.0/
        DATA RAY_L/600D-9/, RAY_W0/0.5D-3/
        DATA PHYS_U0/1.0/,
     &       PHYS_SIGMA/5D-8/, PHYS_SIGMA/5D-14/, PHYS_RHO/1000D0/,
     &       PHYS_CV/4200D0/, PHYS_K/0.61D0/, PHYS_KD/0.61D0/,
     &       PHYS_N0/1.0D0/, PHYS_DNDT/-1.0D0/
     &       TRI_NUM/8000/

        real*8 CFG_VAL
        character*64 CFG_KEY

        open(54,FILE='cell.cfg')
 1101   read(54,*,END=1102) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'L') then
          DIM_L=CFG_VAL
        elseif (CFG_KEY.EQ.'D') then
          DIM_D=CFG_VAL
        elseif (CFG_KEY.EQ.'Lw') then
          DIM_Lw=CFG_VAL
        elseif (CFG_KEY.EQ.'Dw') then
          DIM_Dw=CFG_VAL
        elseif (CFG_KEY.EQ.'Ld') then
          DIM_Ld=CFG_VAL
        elseif (CFG_KEY.EQ.'Dd') then
          DIM_Dd=CFG_VAL
        elseif (CFG_KEY.EQ.'Dr') then
          DIM_Dr=CFG_VAL
        elseif (CFG_KEY.EQ.'Li') then
          DIM_Li=CFG_VAL
        elseif (CFG_KEY.EQ.'Lo') then
          DIM_Lo=CFG_VAL

        elseif (CFG_KEY.EQ.'RAY_L') then
          RAY_L=CFG_VAL
        elseif (CFG_KEY.EQ.'RAY_W0') then
          RAY_W0=CFG_VAL

        elseif (CFG_KEY.EQ.'U0') then
          PHYS_U0=CFG_VAL

        elseif (CFG_KEY.EQ.'SIGMA') then
          PHYS_SIGMA=CFG_VAL
        elseif (CFG_KEY.EQ.'SIGMAd') then
          PHYS_SIGMAD=CFG_VAL
        elseif (CFG_KEY.EQ.'RHO') then
          PHYS_RHO=CFG_VAL
        elseif (CFG_KEY.EQ.'CV') then
          PHYS_CV=CFG_VAL
        elseif (CFG_KEY.EQ.'K') then
          PHYS_K=CFG_VAL
        elseif (CFG_KEY.EQ.'Kd') then
          PHYS_KD=CFG_VAL
        elseif (CFG_KEY.EQ.'N0') then
          PHYS_N0=CFG_VAL
        elseif (CFG_KEY.EQ.'DNDT') then
          PHYS_DNDT=CFG_VAL

        elseif (CFG_KEY.EQ.'TRI_NUM') then
          TRI_NUM=CFG_VAL

        else
          write(*,'(A,A20)')
     *     'warning: unknown parameter in cfg-file: ', CFG_KEY
        endif
        goto 1101
 1102   close(54)
      end
