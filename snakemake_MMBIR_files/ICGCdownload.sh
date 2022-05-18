module unload jdk
module load openjdk/11.0.2

#key=8d6269c8-e16d-4464-8860-96bbec14cb27

~/bin/score-client-2.0.0/bin/score-client --profile $1 download --manifest $2 --output-dir $3

#Before you can actually download controlled access data, you will need to add the ICGC access token in the following file: score-client-x.x.x/conf/application.properties. In the file, you should see a line like shown below. Uncomment it, add your access token, and save the file.
# accessToken=your_collab_access_token
#Assuming that you downloaded a manifest file (manifest.collaboratory.1525977569066.tsv) from ICGC data portal in the previous step, the manifest will contain files from Collaboratory. The following command will download these files to the current directory:





