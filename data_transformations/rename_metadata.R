
#re-name metadata


meta_16S <- meta_16S%>%
  dplyr::rename(O2umol.l = "O.conc..µmol.l.")%>%
  dplyr::rename(NO3_umol.l = "NO3..µmol.l.")%>%
  dplyr::rename(NO2_umol.l = "NO2..µmol.l.")%>%
  dplyr::rename(NH4_umol.l = "NH4..µmol.l..")%>%
  dplyr::rename(PO4_umol.l = "PO4..µmol.l.")%>%
  dplyr::rename(Si_umol.l = "Si..µmol.l.")

meta_18S <- meta_18S%>%
  dplyr::rename(O2umol.l = "O.conc..µmol.l.")%>%
  dplyr::rename(NO3_umol.l = "NO3..µmol.l.")%>%
  dplyr::rename(NO2_umol.l = "NO2..µmol.l.")%>%
  dplyr::rename(NH4_umol.l = "NH4..µmol.l..")%>%
  dplyr::rename(PO4_umol.l = "PO4..µmol.l.")%>%
  dplyr::rename(Si_umol.l = "Si..µmol.l.")
