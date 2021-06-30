# PhlashyNAMe Tutorial

## Introduction 

PhlashyNAMe is a command line tool for downstream analysis Proteomics and Phosphoproteomics data. The tool can be run on linux, Mac and Windows operating systems. What makes PhlashyNAMe unique is it's ability to use phospho-peptide level information to analyse your data. Reactome annotates proteins in specific phosphorylation states as separate entities, therefore these phosphoproteins can take part in different signalling cascades depending on their phosphorylation state. Proteins in these differing modification states are hereby referred to as proteoforms. Taking advantage of this structure, PhlashyNAMe takes an input file of phosphopeptides and maps the phosphorylations onto proteins and complexes according to a set of guidelines explained in detail in the figure below. It works by assigning a confidence score and an abundance score to each mapped protein and complex. 
#### Confidence Score 
This confidence score summarises support for a given phosphorylation state of a protein over its other phosphorylation states based on observations from the data. 

In Reactome, each modified protein is treated as a separate entity to the unmodified protein. In this document these entities are referred to as proteoforms, as they are different versions of the same protein and participate in different parts of the network accordingly. Proteoforms are highlighted in green throughout the  below figures.

![Reactome Proteoforms](https://user-images.githubusercontent.com/9949832/121049429-06383e00-c7fb-11eb-8a4d-e9677ad0a220.png)


Across all proteoforms (green), all recorded phosphorylations are aggregated into one group called Points of Interest (orange) as shown in above. The phosphorylations recorded in the network in this example (PA,PB, and PC) are shown in purple while, phosphorylations found in the data (and not in the network) are shown in red. 

![Support Score](https://user-images.githubusercontent.com/9949832/120878866-98292680-c602-11eb-9e33-aaf8e3549ee1.png)

Looking at the rules developed to score each proteoform, we first score each peptide and take the average of the peptide scores to get the proteoform Confidence score. In rule 1 of the peptide score, a match means at the amino acid in the proteoform, and the peptide have the same status (both phosphorylated or both unphosphorylated). We give a score of 0.5 for a mismatch because that peptide supports the existence of that proteoform, but we know that peptide contains a P.O.I. and will support another proteoform better. Next, looking at Rule 2, we take the average of all possible matches. In Rule 3, we multiply by 0.9 to reflect that the database takes precedent over the data. However, we also reduce the weight of the score to reflect the extra unknown phosphorylation. Finally, in Rule 4, because a peptide mapping perfectly to a modified proteoform is uncommon and of interest we multiply the score by 1.5 to highlight the match. However, we do not highlight this if the proteoform is unmodified as perfect matches would be unmodified peptides (which mostly occur from unspecific binding). 


#### Abundance Score 
Concurrently but separately, abundances (e.g., intensity, SILAC ratio, p-value, log fold change) are mapped onto all proteins and complexes across the network. The process of computing this score is illustrated in the figure below. The default abundance score mapping takes the abundance from the peptide with the highest confidence score (if there are multiple peptides with the same confidence score the average of the abundances are taken). The user can also choose to compute the mean or median abundance for each phosphopeptide mapped to a protein. 

![Abundance score](https://user-images.githubusercontent.com/9949832/120879240-5b126380-c605-11eb-8d3f-6b691f454cfd.png)

When the confidence score is mapped, the abundance score is mapped simultaneously, although the scores remain separate. The user has 4 options for mapping the abundance score as described above. 

## Dependencies 
* Java 8
* R
  * ggplot2 
  * dplyr
  * viridis
  * ggExtra 
  * devtools

## Setting up 

1. Clone this repo (as below) or create a new directory and place the provided scripts 

```
git clone https://github.com/HannahHuckstep/Db_Compare.git
```
2. If you do not have the required depndencies the following `conda` command will create an environment called `PhlashyNAMe` with all the dependencies installed. 

```
conda create --name PhlashyNAMe \
  --channel conda-forge \
  --channel bioconda \
  r=3.6 \
  r-devtools \
  r-ggplot2 \
  r-dplyr \
  r-viridis \
  r-ggExtra \
  openjdk=8
```

3. [OPTIONAL] download the most current version of Reactome and PhosphositePlus. There are currently pre-made integrated and non-integrated databases which can be found in the databases file in this repo for both human and mouse. 

4. Map your data as shown below 

5. Data can then be analysed with useing the number of functions provided below with detailed examples. 

## Tool options and commands 
To start, navigate into the repo directory and type the following command to view all of the tool options: 

```
java -jar jars/ReactoSitePlus.jar -h
```
The first 2 arguments are named -h (or --help) to access the above again (This will be useful later). While the second named option -m (or --mode) is used to specify which function you would like to perform. the current options are listed in the '{}' and then again below in a list form. As an example, the first option is: 
* "CreateDB", takes an  OWL  file  [-iof],  an  output  path  [-op],  an  optional  update  boolean  [-u]  (can  be  T  or  F,  default  is  T),  and the
                         species of graph you'd like to make [-s] (can be human (h) or mouse(m))

So if you would like to create a database you would put the `CreateDB` option after the -m when performing the command. However, this function requires a few parametes to be specified which are shown in the '\[\]' 

For this option you will need to specify:
* the name of the **i**nput **O**WL **f**ile \[-iof\] of [Reactome](https://reactome.org/download/current/biopax.zip) that you would like to build. 
* the **o**ut***p***ut directory where you'd like your graph to reside \[-op\]
* optionally you can specify if you'd like your database to be **u**pdated \[-u\] which can only be set to 'T' (true) or 'F' (false). (updating the database will update secondary UniProt accessions (reffered to here as UniProt ID's) to current UniProt ID's, it will also annotate any deleted or non-human UniProt ID's as deleted) 
* you must also specify which **s**pecies of database you'd like to create \[-s\], which can be set to Human with 'h', 'human', '9606', or Mouse 'm', 'mouse', '10090'

Thus, the final command would look something like this: 
```
java -jar ./path/to/jars/RactoSitePlus.jar -m CreateDB -iof ./path/to/file/Reactome.owl -op ./path/to/graph/ -u T -s h
```

**All comands will start with** ```java -jar ./path/to/jars/RactoSitePlus.jar```

Finally, below all the commands are all of the parameters required to run each command, and an explanation of that they are used for. e.g.,
* --input_owl_file [INPUT_OWL_FILE], -iof [INPUT_OWL_FILE]
                         The OWL file to input

## Mapping Phosphoproteomic Data 

Now we should be ready to map our data! Data can be mapped to pre-built neo4j databases located in the databases folder, or you can build your own up-to-date database as shown above (and again with more detail below). 
This command is shown in the help guide as: 
* "MapPeptides",  takes  an  input  database  [-idb]  and  an  output  path  [-op],   a   file  to  map  onto  the  database  [-idf],  and  the  optional
                         Abundance Score mapping method preferred [-as] ("HighestSupport" is defalut)

Which needs the following parameters specified: 
* the input database directory \[-idb\]
* the output directory for the mapping report \[-op\]
* the file of data you would like to map onto the database \[-idf\]
* The mapping option you'd prefer. There are 5 options 'HighestSupport', 'Max', 'Mean', 'Median', and 'Extreme'
    * Each option is explained above in the abundance score section. The default option is 'HighestSupport'.  



To map the data you can run the following command
``` 
java -jar ./path/to/jars/RactoSitePlus.jar -m MapPeptides -idb ./path/to/graph/ -op ./path/to/output/ -idf ./path/to/data.txt -as HighestSupport
```
                         
In order to map your data you will need an input file with the following 4 columns; 
1. A column where each cell contains a single UniProt ID corresponding to a peptide. *The name of the column cannot contain spaces*. 
    * e.g., leading_razor_protein 
2. A column specifying the modified peptide sequence. *The name of the column cannot contain spaces*. 
    * The modified peptide may be one of 2 different formats. Either \_(ac)AAAITDM(ox)ADLEELSRLS(ph)PLPPGS(ph)PGSAAR\_  or AAAITDMADLEELSRLpSPLPPGpSPGSAAR
3. A column containing the value you would like mapped. *The name of the column cannot contain spaces*. 
    * e.g., the column that holds: p-value, SILAC_ratio, log2_intensity
4. A column you would like the mapped values associated with, such as the experiment name or time. *The name of the column cannot contain spaces*. 
    * e.g., the column that holds:  stimulated, unstimulated, time_point_1

For demonstration purposes we will be mapping phosphoproteomic data generated by [Humphrey et al.](https://www.nature.com/articles/nbt.3327), which can be found in ProteomeXchange under[PXD001792](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001792). The dataset chosen for this tutorial is from the Mouse liver cell line FL38B where 6 replicates were treated with PBS (control) or Insulin. This dataset was processed and select rows taken and made available in the data folder of the corresponding git repository. 

A quick look at the data: 

![SeansData screenshot](https://user-images.githubusercontent.com/9949832/123599616-24f98580-d839-11eb-8784-79f443b9aadd.png)


For this example the parameters we set are: 
1. LRP
2. mod.pep.seq
3. Log2INT
4. Experiment

![Filling out parameters](https://user-images.githubusercontent.com/9949832/123599868-6b4ee480-d839-11eb-81f3-73a37ee65293.png)

The tool will automatically detect the different experiments and map them seperately. In this case, there are 2 experiments 'Control' and 'Insulin'. So after mapping we can look at the mapping report generated for each experiment named in the data named 'PhosphoMappingReport_Control' and 'PhosphoMappingReport_Insulin' to understand how well each experiment mapped. Looking at the Control report: 

![Screen Shot 2021-06-28 at 5 54 39 pm](https://user-images.githubusercontent.com/9949832/123600527-1d86ac00-d83a-11eb-9447-7adbc22e8b08.png)

The first few lines give the statistics from the data, indicating the number of peptides and proteins that mapped etc., while the network statistics indicate the number of things in the database that were mapped to. Following that, the number of proteins and complexes that have been measured in each cellular location are shown in decending order. 

![Screen Shot 2021-06-28 at 5 54 56 pm](https://user-images.githubusercontent.com/9949832/123600545-237c8d00-d83a-11eb-8439-8c4cd5b3b8e3.png)

In addition, the Reactome pathways with the highest proportion of measured proteins and complexes are listed. The proportion measured is shown at the end of each pathway along with the size of the size of the pathway. 



We can also look at our data in relation to qPhos. qPhos is a database holding 554 different experiments accross 137 human cell lines. In order to provide context and understand how well your experimental data mapped, we mapped all 554 experiments to the integrated database and recorded numerous statistics. By running the provided Rscript command( ```Rscript -e "rmarkdown::render('proportionPlots.Rmd')"``` ) in the same directory as the provided 'R' directory you can visualize your experiment in relation to all qPhos experiments as well as other general mapping information. 

![Screen Shot 2021-06-28 at 11 10 53 pm](https://user-images.githubusercontent.com/9949832/123642097-4329aa80-d866-11eb-87d5-3f18f7a91348.png)

This plot depicts the distribution of the proportion of peptides mapped from each of the 554 qPhos experiments and where the Control experiment sits (the red line) in comparison. 

![Screen Shot 2021-06-28 at 11 11 16 pm](https://user-images.githubusercontent.com/9949832/123642081-3e64f680-d866-11eb-92c2-41dea0b25884.png)

This plot depicts the distribution of the proportion of phosphorylation nodes in the database that were mapped to from each of the 554 qPhos experiments and where the Control experiment sits (the red line) in comparison. 

### Visualization 
#### Cytoscape
After the database is mapped to you can write it to a Simple Interaction Format (SIF) file to import into [cytoscape](https://cytoscape.org/) along with an attribute file. The SIF file will be names SIF.sif, which can be renamed but must keep the .sif extension. An example of a SIF file looks like: 

```
21	INPUT	20
23	PHOSPHORYLATION	21
20	OUTPUT	24
```
Which we would read as the node with the id ```21``` is an ```INPUT``` to the reaction node with the id ```20```, or the node with the id ```23``` is a ```PHOSPHORYLATION``` on the node ```21```.

The attribute file will contain all attributes associated with each node. An example of an attribute file looks like: 
```
Node_ID Database_ID Display_Name Type Database_Link Location Status Kinase Transcription_Factor Cell_Surface_Receptor UniProt_Gene_Name Integrated ABUNDANCE_SCORE_wt SUPPORT_SCORE_wt ABUNDANCE_SCORE_stim SUPPORT_SCORE_stim
21 Protein2090 p-S568-MLXIPL Protein http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-HSA-163687.1 nucleoplasm   TRANSCRIPTION_FACTOR  MLXIPL  3.2 0.9 -1.5 0.5
24 SmallMolecule848 Pi SmallMolecule http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-ALL-113550.4 nucleoplasm           
20 BiochemicalReaction1310 [Dephosphorylation of pChREBP (Ser 568) by PP2A] BiochemicalReaction http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=R-HSA-164056.2
23 Protein2090_p-S_568 p-S_568 p_S 568
```

Where the attributes are associated with the node ids found in the SIF file. The first line contains all of the attributes that are found in that database.

Attribute name | Description
---------------|-----------
Node_ID | The unique node number used to link nodes to attributes in cytoscape
Database_ID | The unique node ID given by Reactome 
Display_Name | The name of the node 
Type | The type of node 
Database_Link | The link to the database that has more information on the node. Can link back to Reactome, PhosphoSitePlus or UniProt
Location | The cellular location of the node 
Status | The status of UniProt nodes, can be nothing (if the database was not updated), 'Current', 'Updated', or 'Deleted?' - the question mark indicates further review is needed. The UniProt Id could be old and outdated or not consistent with the species in the database.  
Kinase | If the node is a Kinase this column with have the values 'KINASE' otherwise it will be left blank. These are the lists of kinases for [Human](https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22), and  [Mouse](https://www.uniprot.org/uniprot/?query=keyword:%22Kinase%20[KW-0418]%22&format=list&fil=organism:%22Mus%20musculus%20(Mouse)%20[10090]%22).
Transcription_Factor | If the node is a Transcription Factor this column with have the values 'TRANSCRIPTION_FACTOR' otherwise it will be left blank. These are the lists of Transcription factors for [Human](https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list), and [Mouse](https://www.uniprot.org/uniprot/?query=goa:(%22DNA-binding%20transcription%20factor%20activity%20[3700]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list).
Cell_Surface_Receptor | If the node is a Cell Surface Receptor this column with have the values 'CELL_SURFACE_RECEPTOR' otherwise it will be left blank. These are the lists of Cell surface receptors for [Human](https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22)&format=list) and [Mouse](https://www.uniprot.org/uniprot/?query=goa:(%22cell%20surface%20receptor%20signaling%20pathway%20involved%20in%20cell-cell%20signaling%20[1905114]%22)%20(reviewed:yes%20organism:%22Mus%20musculus%20(Mouse)%20[10090]%22)&format=list).
UniProt_Gene_Name | The gene name associated with that UniProt ID from UniProt.
Integrated | If the node is from PhosphoSitePlus the value will be True. 
ABUNDANCE_SCORE_ | The abundance score mapped from the data for a node. There may be multiple Abundance scores in a single database derived from multiple experiments. Each row labelled ABUNDANCE_SCORE_ will contain the experiment name in the suffix after the last '_'. 
SUPPORT_SCORE_ | The support score mapped from the data for a node. There may be multiple Support scores in a single database derived from multiple experiments. Each row labelled SUPPORT_SCORE_ will contain the experiment name in the suffix after the last '_'. 

To write the SIF and attribute files the java command requires the input database directory and the path to the direcory you'd like the output created in. The command is as follows: 
``` 
java -jar ./path/to/jars/RactoSitePlus.jar -m WriteDBtoSIF -idb ./path/to/graph/ -op ./path/to/output/
```

To then load the SIF file into cytoscape, first open the cytoscape application. Next click the network button (highlighted in red in the figure below).

![Cytoscape Network Button](https://user-images.githubusercontent.com/9949832/120894578-81172280-c65c-11eb-9e3a-d2da0ccb2ff5.png)

Next, navigate to your SIF file and open it. I do not recommend a network view is made at this point as the network is incredibly large and cytoscape often crashes while attempting to make a network view this large. Following this, click the attribute file button (highlighted in red in the figure below). 

![Cytoscape Attribute Button](https://user-images.githubusercontent.com/9949832/120894665-fb47a700-c65c-11eb-8576-4179716395da.png)

Navigate to your attibute file and open it. As in the figure below, make sure the Node_ID column is chosen as the key column. You can also change the value type of the attributes. One change you may like to make is to ensure the ABUNDANCE_SCORE_ and SUPPPORT_SCORE_ columns are set to numeric values. 

![Cytoscape attributes](https://user-images.githubusercontent.com/9949832/120894779-950f5400-c65d-11eb-871a-9d3e8a0c3597.png)

Now your database should be ready to be explored. You may select nodes from the node table (make sure to left click the highlighted nodes and choose 'Select nodes from selected rows') and create a newtork view using the button highlighted in the figure below. 

![Cytoscape new network button](https://user-images.githubusercontent.com/9949832/120895262-8164ed00-c65f-11eb-9c49-383c636fe57f.png)


#### Neo4j
To view and interact with the embedded neo4j databse please follow these instructions: 
1. Download and install the [Neo4j](https://neo4j.com/download-center/#community) (community edition), version 3.5.X as there is a compatibility issue with the later versions.
2. Untar/unzip Neo4j tar/zip file.
3. Copy your current mapped graph directory into a new directory named graph.db. 
    1. You can do this in the command line using this command (for Mac/Linux): ```cp -r path/to/graph/ /graph.db/```
4. Now that you have your graph.db directory, move it into /path/to/neo4j/databases/ and start it with the command ```./path/to/neo4j/bin/neo4j start``` (stop with ``` ./path/to/neo4j/bin/neo4j stop```). It should now be available to view in [localhost](http://localhost:7474/browser/)
    1. Or you can do this proccess manually in the Neo4j desktop app as shown below: 
    2. First open the app and make a new project
    ![Screen Shot 2021-06-06 at 1 42 03 pm](https://user-images.githubusercontent.com/9949832/120911909-582f7580-c6ce-11eb-8b08-114368301918.png)
    4. In the databases section, choose add database. Then, select 'create a local graph'
    ![Screen Shot 2021-06-06 at 1 42 15 pm](https://user-images.githubusercontent.com/9949832/120911922-654c6480-c6ce-11eb-8d66-41ebd101bb5d.png)
    7. Set the name an password to your choosing and create
    ![Screen Shot 2021-06-06 at 1 42 28 pm](https://user-images.githubusercontent.com/9949832/120911946-8614ba00-c6ce-11eb-959c-26bd88788176.png)
    9. Then choose the menu in the top right-hand corner (highlighted in red) and select 'Manage database' 
    ![Screen Shot 2021-06-06 at 1 42 46 pm](https://user-images.githubusercontent.com/9949832/120911966-ba887600-c6ce-11eb-8330-1eaff43221ae.png)
    ![Screen Shot 2021-06-06 at 1 43 02 pm](https://user-images.githubusercontent.com/9949832/120911977-d3912700-c6ce-11eb-8eca-9b06f077016b.png)
    10. Next, select the open terminal button (highlighted in red)
    ![Screen Shot 2021-06-06 at 1 43 14 pm](https://user-images.githubusercontent.com/9949832/120912008-1521d200-c6cf-11eb-9d6c-53db72465164.png)
    12. Once the terminal is open you can copy your graph.db directory into the local directory. The command required is highlighted in red as well as available here (for Mac/Linux) : ```cp -r path/to/graph.db ./```
    ![Screen Shot 2021-06-06 at 1 46 09 pm](https://user-images.githubusercontent.com/9949832/120912068-8c576600-c6cf-11eb-9a8a-47068f9006dd.png)
    13. After the graph.db is in Neo4j, you can press play to start the local database. (highlighted in red) 
    ![Screen Shot 2021-06-06 at 1 47 47 pm](https://user-images.githubusercontent.com/9949832/120912100-c7f23000-c6cf-11eb-984b-06f8094d7806.png)
15. When the database is ready, you can open it to view in your [localhost](http://localhost:7474/browser/)
    ![Screen Shot 2021-06-06 at 1 49 36 pm](https://user-images.githubusercontent.com/9949832/120912133-10115280-c6d0-11eb-83c5-84c648de5d66.png)
16. You can then use the Neo4j query language [Cypher](https://neo4j.com/developer/cypher/). An example of looking a node up by id is below. 
    ![Screen Shot 2021-06-06 at 1 50 58 pm](https://user-images.githubusercontent.com/9949832/120912168-654d6400-c6d0-11eb-952f-4fb634eccb0b.png)
4. if a Windows user,installation and database instrucitions can be found [here](https://neo4j.com/docs/operations-manual/current/installation/windows/) 

*More instructions are available from [Neo4j](https://neo4j.com/docs/operations-manual/current/installation/) and [Reactome](https://reactome.org/dev/graph-database)*


## Analysing Mapped Network 

Now that the network has the data mapped to it, there are a number of ways to analyse the mapped network. 
* Traversal analysis 
* Neighbourhood analysis 
* Find the shortest path between two proteins 
* Minimal Connection Network
* Manual investigation via cytosscape or Neo4j

### Traversal Analysis 

With this function we can look at everything downstream or upstream of a protein of interest. 

TraversalAnalysis, takes in a measured input database [-idb], an output path [-op], a UniProt ID or database ID to look downstream of [-p], the direction of the traversal [-dir], and the experiment name of interest [-en]. Using the example data from earlier these are the function inputs for an analysis looking downstream of the mouse insulin receptor (INSR)(P15208): 

```
java -jar jars/ReactoSitePlus.jar -m TraversalAnalysis -idb ./path/to/graph/ -op ./path/to/output/ -p P15208 -dir downstream -en Control
```

We can look at the resulting report titled "TraversalReport_downstream_P15208.tsv":

![Screen Shot 2021-06-28 at 11 43 54 pm](https://user-images.githubusercontent.com/9949832/123647418-66a32400-d86b-11eb-916b-bd17407c62d0.png)

We see here that this UniProt ID has multiple proteoforms attached to it. Each may be in a different cellular location or have a different profile of modifications. The there a number of statistics reported for the downstream network of each proteoform and are ordered from largest to smallest (in reference to downstream network size). The first line is the variable names that can refer to that particular proteoform, next is the number of nodes found downstream of this proteoform. This statistic includes biochemical reaction nodes, gene nodes etc., as well as nodes that are mappable (proteins and complexes). The remaining statistics are self explanitory. This traversal will traverse all edges excepting small molecule edges, so as to avoid an explosion of things downstream of a single small molecule (such as ATP). 

We can look at the downstream networks in cytoscape by loading the file titled "P15208_downstream.tsv". Once this is loaded into cytoscape and attached to our already loaded measured network, we can filter the network to display the downstream networks of each proteoform. 

![Screen Shot 2021-06-29 at 1 42 07 am](https://user-images.githubusercontent.com/9949832/123665036-42e7da00-d87b-11eb-91eb-ef8cb452c1ea.png)

Looking at the downstream network of INSR (node 55682) in Cytoscape we can start to investigate the differences in expression between experiments directly downstream of our protein of interest. 


This function can also be performed on a particular node ID. This can be any node including reaction or a complex nodes. An example of this function is as follows: 

```
java -jar jars/ReactoSitePlus.jar -m TraversalAnalysis -idb ./path/to/graph/ -op ./path/to/output/ -p 55682 -dir downstream -en Control
```

### Neighbourhood Analysis

With this function we can prioritize neighbourhoods of signalling. In order to prioritize important neighbourhoods, the Empirical False Discovery rate was determined for each mapped neighbourhood. Pre-calculated distributions have been calculated for each neighbourhood of varying sizes in which you can compare your mapped data to. As explained above, qPhos is a database holding 554 different phosphoproteomic experiments accross 137 human cell lines. qPhos was sampled from used to ... ???? bootstrap ? something something.
Can generate your own bkgd distributions with this func : ``` ``` be sure to use the same depth, and approximate experiment size. Also may want to use a super computer as the compute needed is very intensive (get mem stats ~32 GB). 

Using the example data from earlier these are the function inputs: 
```
java -jar jars/ReactoSitePlus.jar -m NeighbourhoodAnalysis
```

We can look at the report:

We can visualize relative statistics in R:

### Shortest Path 

With this function you can look for the shortest path in the network between 2 proteins of interest. 

ShortestPath takes in a measured input database [-idb], an output path [-op], a starting node id [-sid], a ending node id [-eid], and the weight type to be traversed [-ew] (can be either "Abundance" (a) or "Support" (s)). This algorithm will traverse all edge types excepting edges to and from small molecules. Using the example data from earlier we can look at the shortest path between the Insr (P15208) and the Jun protein node seen in the traversal analysis screenshot (node ID 20006):

```
java -jar jars/ReactoSitePlus.jar -m ShortestPath -idb ./path/to/graph/ -op ./path/to/output/ -sid P15208 -eid 20006 -ew a
```
**The function will never be able to find a path between 2 UniProt ID's (or UniProt ID nodes), only between UniProt ID's and node IDs, and node IDs to node IDs** 

We can look at the report titled "ShortestPath_P15208_to_20006.tsv":

![Screen Shot 2021-06-29 at 9 24 34 pm](https://user-images.githubusercontent.com/9949832/123790193-62d1d900-d921-11eb-8914-bfa391ac60bf.png)

This function will report the shortest path between the 2 nodes in every experiment currently mapped on your database. The first line determines if the end node is up- or downstream of the starting node. Next a number of statistics are given. Finally, the path is shown with each node id and the total path weight, followed by the same path but with node names and the weight of each edge beside the reaction nodes.  

We can visualize in cytoscape by loading the file "P15208_to_20006_downstream.tsv	":

![Screen Shot 2021-06-29 at 9 26 36 pm](https://user-images.githubusercontent.com/9949832/123790331-87c64c00-d921-11eb-81e2-3e2f240a3667.png)

### Minimal Connection Network 

With this function you can generate the network of all shortest paths between all nodes (typically reffered to as the shortest-paths network). 

Continuing with the example data the function inputs call would look something like this: 

```
java -jar jars/ReactoSitePlus.jar -m MinimalConnectionNetwork -idb ./path/to/graph/ -op ./path/to/output/ -en Control
```

The report gives general statistics about what is found in the network by looking at "MinimalConnectionNetworkReport_Control.tsv". 

![Screen Shot 2021-06-30 at 7 14 01 pm](https://user-images.githubusercontent.com/9949832/123936119-3dec6d00-d9d8-11eb-88c1-8d0623747645.png)

We can then visualize in cytoscape by loading the file "MinimalConnectionNetworkReport_Control_cytoscape.tsv"

![Screen Shot 2021-06-30 at 7 20 16 pm](https://user-images.githubusercontent.com/9949832/123936327-6a07ee00-d9d8-11eb-8292-0a1d2db1ed32.png)


## Creating a new integrated Database 

To create an embedded integrated database, first download the latest OWL files from [Reactome](https://reactome.org/download/current/biopax.zip) and [PhosphoSitePlus](https://www.phosphosite.org/homeAction). 
1. Download the latest [Reactome](https://reactome.org/download/current/biopax.zip) database.
2. Untar/unzip Neo4j tar/zip file.
    * you can use either Homo_sapiens.owl or Mus_musculus.owl
3. Download the latest [PhosphoSitePlus](https://www.phosphosite.org/homeAction) OWL file. 
    * Navigate to the downloads section, choose 'Datasets from PSP' 
    * Read and agree to terms and conditions 
    * Download 'BioPAX:Kinase-substrate information' 
4. Create the Reactome embedded database as shown above 
    * e.g.,  ```java -jar ./path/to/jars/RactoSitePlus.jar -m CreateDB -iof ./path/to/file/Reactome.owl -op ./path/to/graph/ -u T -s h```
    * Remember to specify wheather or not you'd like the database to be updated and the species 
    * This function downloads many files live from UniProt. Therefore if internet access is unstable this function may crash. You can simply re-enter the command and it will overwrite the current failed output. 
5. Once the Reactome graph is built, you can integrate with PhosphoSitePlus with the 'IntegratePSP' mode. 
    * e.g., ``` java -jar ./path/to/jars/ReactoSitePlus.jar -m IntegratePSP -idb ./path/to/Reactome/Graph/ -op ./path/for/integration/report/output/ -iof ./path/to/PSP.owl```
    * **One key thing to remember is that the orgininal Reactome database that is input into this function is modified and becomes the integrated database.**

## Other Helpful Functions 
* AmountWithLabel
* PrintDatabase
* WriteAllUIDs
* WritePhos
* GetSpecies
* ResetScores
* PrintAllProperties
* qPhosED


