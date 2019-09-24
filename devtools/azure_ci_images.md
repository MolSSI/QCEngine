# QCEngine Azure CI for Private Code 


Engine is a wrapper for many quantum chemistry codes. However, not all of those codes are public, free, or shipped 
in pre-packaged ways (e.g. conda). Testing such codes on public CI suites such as Travis, Circle, or Azure can be a 
security issue as the codes have to be reachable by the CI system. One feature that the Engine developers offer is 
a secure and private Docker repository hosted through 
[Azure Container Registry (ACR)](https://azure.microsoft.com/en-us/services/container-registry/?cdn=disable). This 
allows us to test these codes (with permission) through a mature CI suite that plugs directly into GitHub, without 
exposing the code to the outside world.

The Container Registry is paid for by MolSSI and only the QCArchive developers have access to it.

## Security

The Azure Container Registry (ACR) is private in that only a single developer has access control to the registry. Any 
additional registry managers must be explicitly given access by the owner to take actions, and those users are 
currently limited to a select few QCArchive developers. Although more control can be given, for now, only these 
users have any permission to view, add, or pull images from the private Registry. These accounts are restricted 
to single-user email addresses. The Docker images uploaded to this registry *do not* exist on DockerHub or DockerCloud; 
they *only* live in the private Azure Container Registry.

Only ACR authorized users are allowed to hook the Private Container Registry into the Azure 
Pipelines account. This must be done from the administration panels of Azure DevOps, and even then can only be done 
if said user is first invited (by an administrator) into the DevOps Pipeline Group. Within DevOps, accessing the 
ACR can only be done if the following conditions are met:

* The user with access to the ACR is considered an "Owner" by the ACR security, this permission can only be 
  authorized by the ACR administrator
* The user with access to the ACR is a member of the Azure DevOps Pipeline team which is wanting to queue up the build.
* The user with access to the ACR is an Administrator on the DevOps Pipeline team
* The user with access to the ACR is signed into DevOps and approves using their credentials for the Pipeline team 
  to access the ACR.
* The user with access to the ACR approves specific builds to have access to the ACR. 

Then, and only then, can a Pipeline build access images in the ACR. Further, no unauthorized user can see, manage, 
or change these credentials.

## CI Runtime Security (When the builds trigger)

Engine is, itself, a public code base, even though the codes it wraps may be private. Similarly, the Azure YAML file 
which has the build and test instructions exist inside the Engine code and could be edited by anyone. In order to 
block malicious edits from running and exposing the proprietary code by triggering the CI, the Engine developers only 
allow builds to trigger on code which has been reviewed and approved, and merged into the code base already. No 
outside user can propose a change which will automatically trigger the CI and execute unapproved code.

The trade off to this level of security is that benign PRs which would actually benefit from said testing cannot 
trigger tests on the secure ACR images. The code will be tested on merge into master, and we are working on a better 
manual trigger on request system to let us test PRs which the administrators deem safe. @QCArchiveBot may make an 
appearance as well through a manual trigger phrase.

## Authenticating a new users for ACR and Pipelines

These are instructions for allowing a new user permission to use the private ACR instance paid for by MolSSI/VT. These 
are a one time thing and if you cannot meet these criteria, skip down to section on "Adding a new code to the CI" and 
coordinate with one of the approved users who can upload your custom Docker image (or help you make one) to the private 
ACR.

1. Have an [Azure Portal](https://azure.microsoft.com/en-us/account/) account through a Virginia Tech email address 
   (this is something which has to be authorized by VT so we can allow access to the ACR, this may be relaxed in the 
   future.)
2. Create an [Azure DevOps](https://azure.microsoft.com/en-us/services/devops/?nav=min) account with the same email 
   address. These are different login nodes, even though they both run through Azure. This will be for access to the 
   Pipelines team.
3. Request access to the QCArchive DevOps team. This will require reaching out to one of the QCArchive developers.
4. Request access to the QCArchive ACR. This will require reaching out to one of the QCArchive developers, this is 
   also a smaller group.

*If you only want to upload images to the ACR (i.e. the Pipeline is already connected to the ACR)*

* Request that you be made a "Contributor" in the ACR administration.

*If you want to manage account connections between the ACR and Pipelines*

* Request that you be made an "Owner" in the ACR administration.
* Request that you be made an "Administrator" in the Pipelines team.

## Connecting ACR to Pipelines

See also (its missing details about the security setup)

https://docs.microsoft.com/en-us/azure/devops/pipelines/library/service-endpoints?view=azure-devops&tabs=yaml

*Note: This only has to be done once per ACR (not per image) and Pipelines Team (not per build).*

### Pre-requisites

* Do everything in the previous section
* Be made an "Owner" of ACR and an "Administrator" of the Pipelines team
* Be signed into the account which meets these criteria.

### Linking the services

1. Sign into the account which meets the pre-requisites through the DevOps login.
2. From the DevOps team, go to the group account, and then go to the "Project settings" blade.
    * If you are not an "Administrator" of the Pipelines team, you will not see this button.
3. Go to the "Service Connections" menu.
4. Choose "+ New Service Connection" and then "Docker Registry"
5. Select the "Azure Container Registry" radio button.
6. Choose a "Connection Name" which is easy to recognize and human readable (arbitrary).
    * Make note of this, you will use it when referencing images.
7. Select the "Azure Subscription" drop down, and choose the Azure Subscription which is linked to the ACR you want to 
   make available to this Pipeline group.
    * This will auto-fill if there is only 1 connected Azure Subscription. If nothing appears, the account you are 
    signed in with is missing permissions and likely not associated with the Azure Subscription (e.g. you are using the 
    wrong account).
8. Select the "Azure Container Registry" drop down and follow the authentication before selecting which ACR (if there
   are several) you want to get authorization to. 
9. Click "OK" and this Pipeline will now have access to the ACR through this Service Connection. All authorization 
   will be handled by Azure DevOps invisibly and not exposed in the builds.
    * If the user is not an "Owner" of the ACR, this will throw an authentication error before the window is closed.
    
See also: 
 
https://docs.microsoft.com/en-us/azure/devops/pipelines/library/service-endpoints?view=azure-devops&tabs=yaml#sep-docreg

## Creating and uploading a new image to the ACR

### Pre-requisites

* Docker
* Access to the private software you want to include for CI testing
* (For uploading) Be a "Contributor" (or otherwise Push access) to the ACR in question. See some of the "Authenticating 
  a new users for ACR and Pipelines" for more details.
* (For uploading) Install the [Azure CLI](https://docs.microsoft.com/en-us/cli/azure/install-azure-cli?view=azure-cli-latest)
  
### Creating the Image

1. Start a new [Dockerfile](https://docs.docker.com/engine/reference/builder/) with the following one line header:
    ```Dockerfile
    FROM condaforge/linux-anvil-comp7
    
    ```
2. Build the rest of your docker file as normal to include your code.
3. Use `docker build` but importantly add `-t qcarchive.azurecr.io/{IMAGE NAME}:{TAG}` to the build command.
    * The `{IMAGE NAME}` and `{TAG}` are developer choices, but the repository location is exact. 

See also (for general ACR manipulation). Use the above instructions for making production CI images

https://docs.microsoft.com/en-us/azure/container-registry/container-registry-get-started-portal

### Uploading the image

Note: There is likely better way to authenticate, but this work for now.

See also: 

https://docs.microsoft.com/en-us/azure/container-registry/container-registry-authentication

1. Using the AZ CLI, sign into the QCArchive ACR with `az acr login --name QCArchive` 
    * This targets the repository name from the "Access keys" of the Azure Portal->ACR page inside the container 
    service.
    * You will be prompted with a window to authenticate with your Azure Portal account. Use the one which 
    has access to the ACR.
2. You now have access to ACR to take `docker` based actions against it for 60 minutes by default.
3. Ensure the image you want to push (`docker image ls`) is correctly tagged such that its REPOSITORY location starts 
   with `qcarchive.azurecr.io` and then is given a formal name and TAG as well.
4. `docker push` your image, and it should be pushed to the private ACR.

## Adding a new container to the CI Pipeline

These instructions are for adding new codes to the CI suite.

### Pre-requisites

* ACR and the Pipeline are linked services and you know the `Connection Name` inside the DevOps Pipeline 
  (see "Linking the services").
* The Docker image has been uploaded to the private ACR and you know the image name and tag you want to use 
  (see "Creating and uploading a new image to the ACR")
* You have access to edit the `.azure-pipelines.yml` file of the Engine repository.

### Adding the Custom Image to the YAML Pipeline Container Resource

Assuming you have met the prerequisites, you can request that the CI use the private docker images in its process by 
specifying they [run as a container](https://docs.microsoft.com/en-us/azure/devops/pipelines/process/container-phases?view=azure-devops&tabs=yaml) 
inside the host agent itself. However, we first have to indicate that this container *can* be used by the builds. 
In the `.azure-pipelines.yml` file, there is a header called `resources:` and a sub-header called `containers:` within 
it. This sub-header defines all the containers we want to use and will look like this:

```yaml
resources:
  containers:
    -   container: 'qcengineci' # Arbitrary name for reference in pipelines
        image: qcarchive.azurecr.io/qcengineci:latest # pointer to the Image in our Azure CR
        endpoint: "QCArchive Azure Container Registry" # Name of the Service Connection the pipeline is configured with
    -   container: 'waffleiron' 
        image: qcarchive.azurecr.io/something_about_waffles:maple_syrup
        endpoint: "QCArchive Azure Container Registry" # Name of the Service Connection the pipeline is configured with
```

In this file, you need to specify each custom image in the ACR you want to use as an entry in the list. So for this 
example, we have 2 images, one which we call `qcengineci` and the other `waffleiron`. The names after the 
`container:` directive are arbitrary, but used later, so keep those in mind.

Next we have the `image:` directive which is the pointer to the image in the ACR of the form 
`qcarchive.azurecr.io/{IMAGE_NAME}:{TAG}`. Even if you did not not upload the image yourself, if you know the name, you 
can request the CI fetch it. In these examples, there are two images we reference, both in the `qcarchive.azureacr.io` 
registry, one image called `qcengineci` with the tag `latest`, and one image called `something_about_waffles` with the 
tag `maple_syrup`.

Finally, and most importantly, we have the `endpoint:` directive. This is the human-defined name when the ACR and the 
Pipeline were linked through a "Service connection" in the `Connection Name` field. The specification of this field 
instructs Azure Pipelines to authenticate through that named service when attempting to pull down the image, but only 
in its own Pipeline, i.e. there is no way to reverse engineer that service hook outside of the Pipeline CI, its just a 
string everywhere else. Without this field, Azure Pipeline will try to pull the repository in the `image` 
field, but fail because it has nothing to authenticate against. No passwords/usernames are exposed through this since 
its a string reference to an authentication which can only be seen, created, or edited by administrators though several 
layers of permissions (see "Security").

Once the containers have been requisitioned, they can be used in jobs later.

### Using a custom image in a CI job

To use a previously allocated container (see "Adding the Custom Image to the YAML Pipeline Container Resource"), we 
need to tell the worker Agent to run its job though it. For example:

```yaml
jobs:
- job: proprietary_ci
  displayName: 'CI with Proprietary Code'
  pool:
    vmImage: 'ubuntu-latest'
  container: qcengineci
```

Everything about this should look like a normal Azure Pipelines job, but we have an additional `container` field which 
matches the `container` filed from the `resources: {containers: [{container:...}]` field from above. This effectively 
says that this job, and all actual `steps`, will be run inside this container.

There are some caveats for this. We use Docker images which are derivatives of the Conda-Forge Anvil, which define an 
`ENTRYPOINT` and a `CMD` directive in their `Dockerfile`. Unfortunately, 
[Azure does things which bypasses those directives](https://docs.microsoft.com/en-us/azure/devops/pipelines/process/container-phases?view=azure-devops&tabs=yaml#linux-based-containers), 
so we have to add something in each of our `steps` which tells the container to provision, i.e.

```yaml
  steps:
  - script: |
      source /opt/conda/etc/profile.d/conda.sh
      # The rest of the script...
```

### Running against multiple images

One avenue for future improvements is having multiple images on the ACR we want to reference. For that, we can
[take advantage of the `strategy` directive](https://docs.microsoft.com/en-us/azure/devops/pipelines/process/container-phases?view=azure-devops&tabs=yaml#multiple-jobs) 
of Azure YAML files to auto-queue up multiple jobs.
