QuaVaProt DB docker app

To run the app, first download the data.zip file and extract to the QuaVaProt_db folder, in a folder titled 'data'
Download link: https://figshare.com/s/93e9bcb3fd4c4513c96b

The default app allows searching of both the GDC and COSMIC datasets.
To run the app with COSMIC searching disabled, replace the app.R file in the /shiny-app/ folder with the app.R in the /shiny-app/COSMIC disabled app/ folder

Then build and run the docker image using the command:
docker-compose up --build
