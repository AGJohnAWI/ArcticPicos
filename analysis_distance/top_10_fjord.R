## eukaryotes
#select the top 10species for each fjord

Tanafjord <- meta_18S_nNA.z%>%filter(Fjord == "Tanafjord")
Tanafjord <- Tanafjord$Site

Laksefjord <- meta_18S_nNA.z%>%filter(Fjord == "Laksefjord")
Laksefjord <- Laksefjord$Site

Porsangerfjord <- meta_18S_nNA.z%>%filter(Fjord == "Porsangerfjord")
Porsangerfjord <- Porsangerfjord$Site

Lyngenfjord <- meta_18S_nNA.z%>%filter(Fjord == "Lyngenfjord")
Lyngenfjord <- Lyngenfjord$Site

Balsfjord <- meta_18S_nNA.z%>%filter(Fjord == "Balsfjord")
Balsfjord <- Balsfjord$Site

Lofoten.HE431 <- meta_18S_nNA.z%>%filter(Fjord == "Lofoten.HE431")
Lofoten.HE431 <- Lofoten.HE431$Site

Lofoten.HE33 <- meta_18S_nNA.z%>%filter(Fjord == "Lofoten.HE33")
Lofoten.HE33 <- Lofoten.HE33$Site

Nordfjord <- meta_18S_nNA.z%>%filter(Fjord == "Nordfjord")
Nordfjord <- Nordfjord$Site

Sognefjord <- meta_18S_nNA.z%>%filter(Fjord == "Sognefjord")
Sognefjord <- Sognefjord$Site

Boknafjord <- meta_18S_nNA.z%>%filter(Fjord == "Boknafjord")
Boknafjord <- Boknafjord$Site

OrustT.Fjord <- meta_18S_nNA.z%>%filter(Fjord == "Orust-Tjörn.Fjord")
OrustT.Fjord <- OrustT.Fjord$Site

Van.Mijen.Fjord <- meta_18S_nNA.z%>%filter(Fjord == "Van.Mijen.Fjord")
Van.Mijen.Fjord <- Van.Mijen.Fjord$Site

Isfjord <- meta_18S_nNA.z%>%filter(Fjord == "Isfjord")
Isfjord <- Isfjord$Site

Kongsfjord <- meta_18S_nNA.z%>%filter(Fjord == "Kongsfjord")
Kongsfjord <- Kongsfjord$Site

Woodfjord <- meta_18S_nNA.z%>%filter(Fjord == "Woodfjord")
Woodfjord <- Woodfjord$Site

Wijdefjord <- meta_18S_nNA.z%>%filter(Fjord == "Wijdefjord")
Wijdefjord <- Wijdefjord$Site

Iceland <- meta_18S_nNA.z%>%filter(Fjord == "Iceland")
Iceland <- Iceland$Site

Nordvestfjord <- meta_18S_nNA.z%>%filter(Fjord == "Nordvestfjord.Scoresby.Sund")
Nordvestfjord <- Nordvestfjord$Site

Disco_B <- meta_18S_nNA.z%>%filter(Fjord == "Disco Bay")
Disco_B <- Disco_B$Site



fjord_euk <- euk_r
fjord_euk$Tanafjord <- rowSums(fjord_euk[ , Tanafjord], na.rm = TRUE)
fjord_euk$Laksefjord <- rowSums(fjord_euk[ , Laksefjord], na.rm = TRUE)
fjord_euk$Porsangerfjord <- rowSums(fjord_euk[ , Porsangerfjord], na.rm = TRUE)
fjord_euk$Lyngenfjord <- rowSums(fjord_euk[ , Lyngenfjord], na.rm = TRUE)
fjord_euk$Balsfjord <- rowSums(fjord_euk[ , Balsfjord], na.rm = TRUE)
fjord_euk$Lofoten.HE431 <- rowSums(fjord_euk[ , Lofoten.HE431], na.rm = TRUE)
fjord_euk$Lofoten.HE33 <- rowSums(fjord_euk[ , Lofoten.HE33], na.rm = TRUE)
fjord_euk$Nordfjord <- rowSums(fjord_euk[ , Nordfjord], na.rm = TRUE)
fjord_euk$Sognefjord <- rowSums(fjord_euk[ , Sognefjord], na.rm = TRUE)
fjord_euk$Boknafjord <- rowSums(fjord_euk[ , Boknafjord], na.rm = TRUE)
fjord_euk$OrustT.Fjord <- rowSums(fjord_euk[ , OrustT.Fjord], na.rm = TRUE)
fjord_euk$Van.Mijen.Fjord <- rowSums(fjord_euk[ , Van.Mijen.Fjord], na.rm = TRUE)
fjord_euk$Isfjord <- rowSums(fjord_euk[ , Isfjord], na.rm = TRUE)
fjord_euk$Kongsfjord <- rowSums(fjord_euk[ , Kongsfjord], na.rm = TRUE)
fjord_euk$Woodfjord <- rowSums(fjord_euk[ , Woodfjord], na.rm = TRUE)
fjord_euk$Wijdefjord <- rowSums(fjord_euk[ , Wijdefjord], na.rm = TRUE)
fjord_euk$Iceland <- rowSums(fjord_euk[ , Iceland], na.rm = TRUE)
fjord_euk$Nordvestfjord <- rowSums(fjord_euk[ , Nordvestfjord], na.rm = TRUE)
fjord_euk$Disco_B <- rowSums(fjord_euk[ , Disco_B], na.rm = TRUE)

euk_fjord_total <- fjord_euk 
euk_fjord_clr <- clr(fjord_euk,1)
euk_fjord_clr_sums <- euk_fjord_clr
euk_fjord_clr_re <-  euk_fjord_clr%>%dplyr::slice_max(order_by = Porsangerfjord, n = 10)
#change to the fjord of interest

#euk_No_glacier_top10 <- taxonomy_18S%>%dplyr::filter(rownames(taxonomy_18S) %in% rownames(euk_fjord_clr_re))


#prokaryotes

#select the top 10species for each fjord

Tanafjord <- meta_16S_nNA.z%>%filter(Fjord == "Tanafjord")
Tanafjord <- Tanafjord$Site

Laksefjord <- meta_16S_nNA.z%>%filter(Fjord == "Laksefjord")
Laksefjord <- Laksefjord$Site

Porsangerfjord <- meta_16S_nNA.z%>%filter(Fjord == "Porsangerfjord")
Porsangerfjord <- Porsangerfjord$Site

Lyngenfjord <- meta_16S_nNA.z%>%filter(Fjord == "Lyngenfjord")
Lyngenfjord <- Lyngenfjord$Site

Balsfjord <- meta_16S_nNA.z%>%filter(Fjord == "Balsfjord")
Balsfjord <- Balsfjord$Site

Lofoten.HE431 <- meta_16S_nNA.z%>%filter(Fjord == "Lofoten.HE431")
Lofoten.HE431 <- Lofoten.HE431$Site

Lofoten.HE33 <- meta_16S_nNA.z%>%filter(Fjord == "Lofoten.HE33")
Lofoten.HE33 <- Lofoten.HE33$Site

Nordfjord <- meta_16S_nNA.z%>%filter(Fjord == "Nordfjord")
Nordfjord <- Nordfjord$Site

Sognefjord <- meta_16S_nNA.z%>%filter(Fjord == "Sognefjord")
Sognefjord <- Sognefjord$Site

Boknafjord <- meta_16S_nNA.z%>%filter(Fjord == "Boknafjord")
Boknafjord <- Boknafjord$Site

OrustT.Fjord <- meta_16S_nNA.z%>%filter(Fjord == "Orust-Tjörn.Fjord")
OrustT.Fjord <- OrustT.Fjord$Site

Van.Mijen.Fjord <- meta_16S_nNA.z%>%filter(Fjord == "Van.Mijen.Fjord")
Van.Mijen.Fjord <- Van.Mijen.Fjord$Site

Isfjord <- meta_16S_nNA.z%>%filter(Fjord == "Isfjord")
Isfjord <- Isfjord$Site

Kongsfjord <- meta_16S_nNA.z%>%filter(Fjord == "Kongsfjord")
Kongsfjord <- Kongsfjord$Site

Woodfjord <- meta_16S_nNA.z%>%filter(Fjord == "Woodfjord")
Woodfjord <- Woodfjord$Site

Wijdefjord <- meta_16S_nNA.z%>%filter(Fjord == "Wijdefjord")
Wijdefjord <- Wijdefjord$Site

Iceland <- meta_16S_nNA.z%>%filter(Fjord == "Iceland")
Iceland <- Iceland$Site

Nordvestfjord <- meta_16S_nNA.z%>%filter(Fjord == "Nordvestfjord.Scoresby.Sund")
Nordvestfjord <- Nordvestfjord$Site

Disco_B <- meta_16S_nNA.z%>%filter(Fjord == "Disco Bay")
Disco_B <- Disco_B$Site



fjord_prok <- prok_r
fjord_prok$Tanafjord <- rowSums(fjord_prok[ , Tanafjord], na.rm = TRUE)
fjord_prok$Laksefjord <- rowSums(fjord_prok[ , Laksefjord], na.rm = TRUE)
fjord_prok$Porsangerfjord <- rowSums(fjord_prok[ , Porsangerfjord], na.rm = TRUE)
fjord_prok$Lyngenfjord <- rowSums(fjord_prok[ , Lyngenfjord], na.rm = TRUE)
fjord_prok$Balsfjord <- rowSums(fjord_prok[ , Balsfjord], na.rm = TRUE)
fjord_prok$Lofoten.HE431 <- rowSums(fjord_prok[ , Lofoten.HE431], na.rm = TRUE)
fjord_prok$Lofoten.HE33 <- rowSums(fjord_prok[ , Lofoten.HE33], na.rm = TRUE)
fjord_prok$Nordfjord <- rowSums(fjord_prok[ , Nordfjord], na.rm = TRUE)
fjord_prok$Sognefjord <- rowSums(fjord_prok[ , Sognefjord], na.rm = TRUE)
fjord_prok$Boknafjord <- rowSums(fjord_prok[ , Boknafjord], na.rm = TRUE)
fjord_prok$OrustT.Fjord <- rowSums(fjord_prok[ , OrustT.Fjord], na.rm = TRUE)
fjord_prok$Van.Mijen.Fjord <- rowSums(fjord_prok[ , Van.Mijen.Fjord], na.rm = TRUE)
fjord_prok$Isfjord <- rowSums(fjord_prok[ , Isfjord], na.rm = TRUE)
fjord_prok$Kongsfjord <- rowSums(fjord_prok[ , Kongsfjord], na.rm = TRUE)
fjord_prok$Woodfjord <- rowSums(fjord_prok[ , Woodfjord], na.rm = TRUE)
fjord_prok$Wijdefjord <- rowSums(fjord_prok[ , Wijdefjord], na.rm = TRUE)
fjord_prok$Iceland <- rowSums(fjord_prok[ , Iceland], na.rm = TRUE)
fjord_prok$Nordvestfjord <- rowSums(fjord_prok[ , Nordvestfjord], na.rm = TRUE)
fjord_prok$Disco_B <- rowSums(fjord_prok[ , Disco_B], na.rm = TRUE)

prok_fjord_total <- fjord_prok
prok_fjord_clr <- clr(fjord_prok,1)

prok_fjord_clr_re <-  prok_fjord_clr%>%dplyr::slice_max(order_by = Nordvestfjord, n = 10)
#change to the fjord of interest

#p_top10 <- taxonomy_16S%>%dplyr::filter(rownames(taxonomy_16S) %in% rownames(prok_fjord_clr_re))


