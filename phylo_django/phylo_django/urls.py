"""URL configuration for phylo_django project."""
from django.urls import include, path

urlpatterns = [
    path('', include('phylo_app.urls')),
]
