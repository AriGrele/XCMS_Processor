source('../code/packages.R')
source('../code/functions.R')

gnps_features        = fread('gnps_output/nf_output/clustering/featuretable_reformated.csv')|>
  transform(pcgroup  = fread('xcms_output/pseudospectra.csv')$pcgroup,
            ID       = `row ID`,
            `#Scan#` = fread('gnps_output/nf_output/networking/clustersummary_with_network.tsv')$`#Scan#`)|>
  merge(fread('gnps_output/nf_output/library/merged_results_with_gnps.tsv'),by='#Scan#',all=T)|>
  as.data.frame()

gnps_network = fread('gnps_output/nf_output/networking/merged_pairs.tsv')|>
  as.data.frame()

g = graph_from_data_frame(gnps_network, directed = F)
gnps_network$cluster = components(g)$membership[match(gnps_network$CLUSTERID1, V(g)$name)]

gnps_features = data.frame('ID'  = c(gnps_network[,1],gnps_network[,2]),
                           'cluster' = rep(gnps_network$cluster,2))|>
  unique()|>
  merge(gnps_features, by = 'ID',all=T)|>
  within({
    cluster[is.na(cluster)] = 0
    pcgroup = paste(cluster,pcgroup,sep='-')})

gnps_features = cbind(gnps_features[,c('ID','cluster','pcgroup','Compound_Name','superclass','class','subclass','MQScore','row retention time','row m/z','SharedPeaks','MZErrorPPM')],
                      gnps_features[grepl('mzML',names(gnps_features))])|>
  within({
    check = !((SharedPeaks>=4)|(MQScore>.9)|(MZErrorPPM<=30)) #remove annotations for low confidence compounds
    Compound_Name[check] = NA
    superclass[check]    = NA
    class[check]         = NA
    subclass[check]      = NA
    MQScore[check]       = NA
    SharedPeaks[check]   = NA
    MZErrorPPM[check]    = NA})

pcgroups = unique(gnps_features$pcgroup)
for(group in pcgroups){
  bar(match(group,pcgroups)/length(pcgroups))
  
  use   = gnps_features$pcgroup == group
  sub   = gnps_features[use,]
  check = which(sub$MQScore == max(sub$MQScore,na.rm=T)) 
  if(length(check>0)){
    gnps_features$Compound_Name[use] = gnps_features$Compound_Name[use][check] #terrible, but gets the job done
    gnps_features$superclass[use]    = gnps_features$superclass[use][check]
    gnps_features$class[use]         = gnps_features$class[use][check]
    gnps_features$subclass[use]      = gnps_features$subclass[use][check]
    gnps_features$MQScore[use]       = gnps_features$MQScore[use][check]
    gnps_features$SharedPeaks[use]   = gnps_features$SharedPeaks[use][check]
    gnps_features$MZErrorPPM[use]    = gnps_features$MZErrorPPM[use][check]}}

gnps_features = replace(gnps_features,gnps_features=='N/A',NA)
gnps_features[c('Compound_Name','superclass','class','subclass')] = 
  replace(gnps_features[c('Compound_Name','superclass','class','subclass')],is.na(gnps_features[c('Compound_Name','superclass','class','subclass')]),'No annotation')

write.csv(gnps_features, 'formatted chemistry/peaks.csv',row.names=F)
####network plot data####
node_attrs = (subset(gnps_features, ID%in%unique(c(gnps_network[,1],gnps_network[,2])))|>
                transform(id = ID))[c('id','superclass')]|>unique()

write.csv(gnps_network, 'formatted chemistry/network_edges.csv',row.names=F)
write.csv(node_attrs,   'formatted chemistry/network_nodes.csv',row.names=F)

g = graph_from_data_frame(gnps_network, directed = T, vertices = node_attrs)
(ggraph(g, layout = 'stress')+ 
    geom_edge_link()+
    geom_node_point(aes(color=(superclass)),size=2))|>
  print()
