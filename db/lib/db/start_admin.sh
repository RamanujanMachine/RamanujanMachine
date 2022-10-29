docker rm docker_admin
docker pull dpage/pgadmin4
docker run --name docker_admin -p 80:80 -e "PGADMIN_DEFAULT_EMAIL=amir@test.com" -e "PGADMIN_DEFAULT_PASSWORD=123456" -d dpage/pgadmin4
