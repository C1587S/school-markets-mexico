# Schooling Markets in Mexico

In Schooling Markets (SM), education services are exchanged. The supply side is composed of school institutions competing to enroll students, while the demand side is formed by students' parents who choose where to enroll their children. As Mizala and Urquiola (2013) highlighted, a critical factor that likely determines the impact of competition within markets is explained by whether parents are informed about schools' effectiveness or value-added and how it influences their school choices. Public policy programs could enhance schools' competition over education quality by increasing parents' information on this field.

However, SM identification is not a trivial task and could include different definitions and approaches. This work proposes a novel methodology for detecting school markets at elementary and high school levels using geospatial analysis for calculating Commuting Zones (CZ), and graph algorithms for community detection using student migrations data in Mexico between 2006 and 2013. 

This repository includes the necessary files to reproduce our results. A more detailed description on the methodology can be found at this [link](https://c1587s.github.io/scadavidsanchez/#/research/school_markets).


### Folder Structure

 ```bash
school-markets-mexico/
├── data
│   ├── buffers
│   │   ├── df_buffers_10kms_prim.rds
│   │   ├── df_buffers_10kms_sec.rds
│   │   ├── df_buffers_5kms_prim.rds
│   │   └── df_buffers_5kms_sec.rds
│   ├── migration
│   │   ├── agregado_dist_prim.rds
│   │   └── agregado_dist_sec.rds
│   └── schools
│       ├── db_master.rds
│       └── master_db_filtered.rds
├── README.md
└── scripts
    ├── App.R
    ├── communities_high_school.R
    ├── communities_primary_school.R
    ├── commuting_zones_high_school.R
    ├── commuting_zones_primary_school.R
    └── utils.R

5 directories, 15 files
 ``` 


### References

Mizala, A., & Urquiola, M. (2013). School markets: The impact of information approximating schools' effectiveness. Journal of Development Economics, 103, 313-335.