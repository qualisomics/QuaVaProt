QuaVaProt DB MutaQuant edition docker app

To run the app, first download the data.zip file and extract to the QuaVaProt_mutaquant folder, in a folder titled 'data' 
Download link: https://figshare.com/s/56569878c75622e5f849

The folder structure should be as so:

-Quavaprot_mutaquant/
├── -data/
│   ├── -Gapp_database_1.csv
│   ├── -Gapp_database_2.csv
│   └── -Gapp_database_3.csv
├── -shiny-app/
│   └── -app.R
├── -Dockerfile
├── -docker-compose.yml
└── -nginx.conf

Then build and run the docker image using the command: docker-compose up --build
