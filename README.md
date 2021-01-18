# README #

Markov chain based model of pedestrian behaviour from: Shelat, S. (2017). Developing an Integrated Pedestrian Behaviour Model for Office Buildings. (MSc), Delft University of Technology. 

Simulating strategic, tactical, and operational level pedestrian behaviour. This MSc disseration was conducted in the context of connected lighting in office environments.

If you find this useful for your work, please cite the thesis or the associated paper (http://dx.doi.org/10.17815/CD.2020.78). See abstract for the tactical level model below.

### How to: ###

* Building layouts described in the thesis (appendix) have to be created in a drafting software to create a .shp file. See getFiles.m
* Change simulation settings for all behaviour levels in getSettings.m
* Run main file.

### Abstract from conference paper (only tactical level) ###

As the number of people working in office buildings increases, there is an urgent need to improve building services, such as lighting and temperature control, within these buildings to increase energy efficiency and well-being of occupants. A pedestrian behaviour model that simulates office occupants’ movements and locations can provide the high spatial and temporal resolution data required for the testing, evaluation, and optimization of these control systems. However, since most studies in pedestrian research focus on modelling specific actions at the operational level or target situations where movement schedules do not have to modelled, a pedestrian behaviour model that can simulate complex situations over long time periods is missing. Therefore, this paper proposes a tactical level model to generate occupant movement patterns in office buildings. The Markov-chain activity-based model proposed here is data parsimonious, flexible in accepting different levels of information, and can produce high resolution output. The mathematical properties of the methodology are analyzed to understand their impact on the final results. Finally, the tactical level pedestrian behaviour model is face validated using a case study of an imaginary office with a simple layout.