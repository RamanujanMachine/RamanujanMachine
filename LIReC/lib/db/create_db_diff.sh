docker exec -i ramanujan_db psql -U postgres postgres < create_db_diff.sql


#sqlacodegen postgresql://postgres:123456@localhost:5432/ramanujan --outfile models.py
