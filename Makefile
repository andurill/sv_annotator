CUR_DIR=`pwd`
run-dev:
	. variant_lookup/bin/activate && \
	export PYTHONPATH=/home/jayakumg/software/lib/gl_variant_lookup/dmp-commons/:/home/jayakumg/software/lib/gl_variant_lookup/flask/:/home/jayakumg/software/lib/gl_variant_lookup/lib
	cd $(CUR_DIR)/server && \
	gunicorn --bind=norma.mskcc.org:8087 --workers=2 --log-file error.log --reload app:dev_app

virtual-environment:
	/dmp/resources/prod/tools/system/python/production/bin/virtualenv variant_lookup
	. variant_lookup/bin/activate 

deploy-dev:
	. variant_lookup/bin/activate && \
	pip install -r requirements.txt
	export PYTHONPATH=/home/jayakumg/software/lib/gl_variant_lookup/dmp-commons/:/home/jayakumg/software/lib/gl_variant_lookup/flask/:/home/jayakumg/software/lib/gl_variant_lookup/lib
	cd $(CUR_DIR)/server && \
	gunicorn --bind=norma.mskcc.org:8087 --workers=2 --log-file error.log --reload app:dev_app
