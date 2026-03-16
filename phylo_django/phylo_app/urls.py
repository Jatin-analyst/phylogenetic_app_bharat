from django.urls import path
from . import views

urlpatterns = [
    path('', views.index_view, name='index'),
    path('upload/', views.upload_view, name='upload'),
    path('generate/', views.generate_view, name='generate'),
    path('download/<str:format>/', views.download_view, name='download'),
]
