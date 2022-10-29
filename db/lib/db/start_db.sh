docker rm ramanujan_db
docker pull postgres
docker run --name ramanujan_db -v /home/amir/projects/current_projects/ramanujanpriv/DB/data:/var/lib/postgresql/data -e POSTGRES_PASSWORD=123456 -p 5432:5432 -d postgres
