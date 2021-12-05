
require(BBmisc)
info <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/info_table_drifter.csv", sep = ";")

#sort according to the number
info.sort <- with(info, info[order(info$no),])

drifter_1_month <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/a4_fwd_cora_cm_day_31.csv", head = F)
drifter_3_month <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/a4_fwd_cora_cm_day_91.csv", head = F)
drifter_6_month <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/a4_fwd_cora_cm_day_181.csv", head = F)
drifter_1year <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/a4_fwd_cora_cm_day_361.csv", head = F)
drifter_5year <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/a4_tracmass2/a4_fwd_cora_cm_day_1821.csv", head = F)

rownames(drifter_1_month) <- info.sort$Station
colnames(drifter_1_month) <- info.sort$Station

#drifter_1_month[drifter_1_month == -999.0000] <- 0
drifter_1_month_p <- drifter_1_month#1/drifter_1_month[] *(-1)

drifter_3_month_p <- drifter_3_month #1/drifter_3_month[] *(-1)

#drifter_6_month[drifter_6_month == -999.0000] <- 0
drifter_6_month_p <- drifter_6_month  #1/drifter_6_month[] *(-1)

drifter_1year_p <-  drifter_1year#1/drifter_1year[] *(-1)

drifter_5year_p <- drifter_5year#1/drifter_5year[] *(-1)

#4proks

drifter_1m_prok <- drifter_1_month_p
rownames(drifter_1m_prok) <- info.sort$Site_prok
colnames(drifter_1m_prok) <- info.sort$Site_prok

drifter_3m_prok <- drifter_3_month_p
rownames(drifter_3m_prok) <- info.sort$Site_prok
colnames(drifter_3m_prok) <- info.sort$Site_prok

drifter_6m_prok <- drifter_6_month_p
rownames(drifter_6m_prok) <- info.sort$Site_prok
colnames(drifter_6m_prok) <- info.sort$Site_prok

drifter_1y_prok <- drifter_1year_p
rownames(drifter_1y_prok) <- info.sort$Site_prok
colnames(drifter_1y_prok) <- info.sort$Site_prok

drifter_5y_prok <- drifter_5year_p
rownames(drifter_5y_prok) <- info.sort$Site_prok
colnames(drifter_5y_prok) <- info.sort$Site_prok

#4euk

drifter_1m_euk <- drifter_1_month_p
rownames(drifter_1m_euk) <- info.sort$Site_euk
colnames(drifter_1m_euk) <- info.sort$Site_euk

drifter_3m_euk <- drifter_3_month_p
rownames(drifter_3m_euk) <- info.sort$Site_euk
colnames(drifter_3m_euk) <- info.sort$Site_euk

drifter_6m_euk <- drifter_6_month_p
rownames(drifter_6m_euk) <- info.sort$Site_euk
colnames(drifter_6m_euk) <- info.sort$Site_euk

drifter_1y_euk <- drifter_1year_p
rownames(drifter_1y_euk) <- info.sort$Site_euk
colnames(drifter_1y_euk) <- info.sort$Site_euk

drifter_5y_euk <- drifter_5year_p
rownames(drifter_5y_euk) <- info.sort$Site_euk
colnames(drifter_5y_euk) <- info.sort$Site_euk


