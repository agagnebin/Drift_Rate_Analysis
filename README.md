Anna Gagnebin (Sacramento State Physics, 2023)

Mentor: Dr. David DeBoer

Berkeley SETI Research Center, Summer 2021

The goal of this project was to develop a Jupyter Notebook to help a user visually analyze the drift of a potential signal of interest to verify whether it: 
- may be extraterrestrial in origin
- could have originated from a source located in the direction in which the signal was observed
in order to simplify the search for extraterrestrial intelligence.

The outline of this repository:
1. Code
- A Jupyter Notebook containing:
  - code that can be easily used to visualize the drift of a potential signal from any source on the sky observed at any location on the Earth
  - a detailed walkthrough of the code
  - an overview of the background research I did this summer on coordinate systems and how to use astropy for drift rate analysis
- Files written by Dave DeBoer that are needed to run the Notebook
2. Presentations
- Research poster from the 2021 Assembly of the Order of the Octopus
- Presentation from the Cal-Bridge Summer Research Symposium

There are a number of steps which may be taken following the end of this summer:
- Next steps:
  - In the Notebook:
   - Edit the Notebook so it may be run independently of the two python files written by Dr. DeBoer
   - Investigate a way to connect to JPL Horizons through the Notebook instead of imbedding a link to the website
     - I tried to use telnetlib, but I couldn't get it to work. The next step would be to look into astroquery
   - Edit the way noise is added into the waterfall plots so they would look more like a real signal
  - In general:
   - Edit the code written by Dr. DeBoer that deals with the other potential causes of drift (outlined in the Notebook)
      - Possibly create a new Notebook for each

<i> Funding for this project was provided by the National Science Foundation and Breakthrough Listen </i>
